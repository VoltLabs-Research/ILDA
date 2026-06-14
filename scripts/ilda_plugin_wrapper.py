from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BINARY_NAME = "ilda"
PLUGIN_REPO_DIRNAME = "ilda"
ENV_BINARY_OVERRIDE = "VOLT_ILDA_BINARY"
REQUIRED_OUTPUTS = ["_bonds.parquet"]
LOG_TAG = "ilda-plugin"

SCRIPT_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = Path(os.environ.get("PLUGIN_PROJECT_DIR", SCRIPT_DIR.parent)).resolve()
PLUGINS_ROOT = PLUGIN_ROOT.parent if PLUGIN_ROOT.parent.name == "plugins" else None
EMBEDDED_LOADER = (PLUGIN_ROOT / "lib/ld-linux-x86-64.so.2").resolve()
EMBEDDED_LIBRARY_DIR = (PLUGIN_ROOT / "lib").resolve()


class WrapperError(RuntimeError):
    pass


def _resolve_binary(binary_name: str, env_var: str, repo_subdir: str) -> Path:
    env_value = os.environ.get(env_var, "").strip()
    candidates: list[Path] = []
    if env_value:
        candidates.append(Path(env_value))
    candidates.append(PLUGIN_ROOT / "bin" / binary_name)
    if PLUGINS_ROOT is not None:
        candidates.extend([
            PLUGINS_ROOT / repo_subdir / "build/build/Release" / binary_name,
            PLUGINS_ROOT / repo_subdir / "build/Release" / binary_name,
            PLUGINS_ROOT / repo_subdir / "build-local/build/Release" / binary_name,
            PLUGINS_ROOT / repo_subdir / "build-manual/build/Release" / binary_name,
        ])
    which = shutil.which(binary_name)
    if which:
        candidates.append(Path(which))
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()
    probed = "\n".join(f"  - {c}" for c in candidates)
    raise WrapperError(f"No pude resolver binario {binary_name}. Paths probados:\n{probed}")


def _resolve_embedded_runtime_command(command: list[str]) -> list[str]:
    if not command:
        return command
    if not EMBEDDED_LOADER.exists() or not EMBEDDED_LIBRARY_DIR.exists():
        return command
    binary_path = Path(command[0]).resolve()
    try:
        binary_path.relative_to(PLUGIN_ROOT)
    except ValueError:
        return command
    if binary_path == EMBEDDED_LOADER:
        return command
    return [
        str(EMBEDDED_LOADER),
        "--library-path", str(EMBEDDED_LIBRARY_DIR),
        str(binary_path),
        *command[1:],
    ]


def _run(command: list[str]) -> None:
    command = _resolve_embedded_runtime_command(command)
    sys.stderr.write(f"[{LOG_TAG}] {' '.join(shlex.quote(part) for part in command)}\n")
    sys.stderr.flush()
    completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=None)
    if completed.returncode != 0:
        raise WrapperError(f"El comando fallo con exit code {completed.returncode}: {command[0]}")


def _require_outputs(output_base: str) -> None:
    for suffix in REQUIRED_OUTPUTS:
        expected = Path(f"{output_base}{suffix}")
        if not expected.exists():
            raise WrapperError(f"Falta el archivo requerido: {expected}")


_VOLT_RUNTIME_FLAGS_WITH_VALUE = {
    "--selectedTimesteps",
    "--selected-timesteps",
}


def _filter_runtime_flags(args: list[str]) -> list[str]:
    filtered: list[str] = []
    i = 0
    while i < len(args):
        token = str(args[i])
        if token in _VOLT_RUNTIME_FLAGS_WITH_VALUE:
            i += 2
            continue
        filtered.append(token)
        i += 1
    return filtered


def _ensure_executable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        if not (mode & 0o111):
            path.chmod(mode | 0o755)
    except OSError:
        pass


def _has_flag(args: list[str], flag: str) -> bool:
    return any(a == flag for a in args)


def _append_flag_if_missing(args: list[str], flag: str, value: str | None) -> list[str]:
    if not value or _has_flag(args, flag):
        return args
    return [*args, flag, value]


def _artifact_path(artifacts: dict, name: str) -> str | None:
    artifact = artifacts.get(name)
    if not isinstance(artifact, dict):
        return None
    path = artifact.get("path")
    return str(path) if path else None


def _find_parent_grain_artifacts(args: list[str]) -> tuple[str, str] | None:
    """Standalone runs: derive upstream artifacts from the output base, the way
    grain-segmentation writes them when chained with outputPathMode:parent."""
    if len(args) < 2:
        return None
    output_base = args[1]
    grain_atoms = f"{output_base}_atoms.parquet"
    grains = f"{output_base}_grains.parquet"
    if all(Path(path).exists() for path in (grain_atoms, grains)):
        return grain_atoms, grains
    return None


def _ensure_grain_artifacts(args: list[str]) -> list[str]:
    if _has_flag(args, "--grain-atoms") and _has_flag(args, "--grains"):
        return args
    parent_artifacts = _find_parent_grain_artifacts(args)
    if parent_artifacts is not None:
        grain_atoms, grains = parent_artifacts
        args = _append_flag_if_missing(args, "--grain-atoms", grain_atoms)
        args = _append_flag_if_missing(args, "--grains", grains)
    return args


def _apply_pipeline_input_overrides(args: list, config: dict) -> list:
    artifacts = config.get("pluginInputArtifacts") or config.get("plugin_input_artifacts")
    if not isinstance(artifacts, dict):
        return args
    next_args = list(args)
    next_args = _append_flag_if_missing(next_args, "--grain-atoms", _artifact_path(artifacts, "grainAtoms"))
    next_args = _append_flag_if_missing(next_args, "--grains", _artifact_path(artifacts, "grains"))
    return next_args


def _run_binary_with_args(args: list[str]) -> dict:
    args = _filter_runtime_flags([str(a) for a in args])
    if len(args) < 2:
        raise WrapperError("Se esperaban al menos 2 argumentos: <input_dump> <output_base>")
    args = _ensure_grain_artifacts(args)

    output_base = args[1]
    Path(output_base).parent.mkdir(parents=True, exist_ok=True)
    binary = _resolve_binary(BINARY_NAME, ENV_BINARY_OVERRIDE, PLUGIN_REPO_DIRNAME)
    _ensure_executable(binary)
    if EMBEDDED_LOADER.exists():
        _ensure_executable(EMBEDDED_LOADER)
    command = [str(binary), *args]
    _run(command)
    _require_outputs(output_base)
    return {
        "ok": True,
        "outputBase": output_base,
        "binary": str(binary),
        "outputs": [f"{output_base}{suffix}" for suffix in REQUIRED_OUTPUTS],
    }


def process(frame, config):
    del frame
    if not isinstance(config, dict):
        raise WrapperError("config debe ser un dict")
    args = config.get("args")
    if not isinstance(args, list):
        raise WrapperError("config['args'] debe ser una lista de strings")
    args = _apply_pipeline_input_overrides(args, config)
    return _run_binary_with_args(args)


def _main_cli() -> int:
    _run_binary_with_args(sys.argv[1:])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(_main_cli())
    except WrapperError as error:
        sys.stderr.write(f"[{LOG_TAG}] error: {error}\n")
        raise SystemExit(1)
