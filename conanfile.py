from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class IldaConan(ConanFile):
    name = "ilda"
    version = "1.0.2"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    requires = (
        "opendxa/[>=2.0]",
        "coretoolkit/[>=2.0]",
        "structure-identification/[>=2.0]",
        "common-neighbor-analysis/[>=2.0]",
        "polyhedral-template-matching/[>=2.0]",
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
        "duckdb/1.4.3",
    )
    default_options = {
        "duckdb/*:shared": False,
        "duckdb/*:with_parquet": True,
        "duckdb/*:with_httpfs": False,
        "duckdb/*:with_json": False,
        "duckdb/*:with_icu": False,
        "duckdb/*:with_tpch": False,
        "duckdb/*:with_tpcds": False,
        "duckdb/*:with_fts": False,
        "duckdb/*:with_inet": False,
        "duckdb/*:with_excel": False,
        "duckdb/*:with_autocomplete": False,
        "duckdb/*:with_visualizer": False,
        "duckdb/*:with_sqlsmith": False,
        "duckdb/*:with_shell": False,
    }
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "scripts/*", "plugin.json", "README.md", ".gitignore"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_target_name", "ilda::ilda")
        self.cpp_info.includedirs = ["include"]
        self.cpp_info.libdirs = ["lib"]
        self.cpp_info.bindirs = ["bin"]
        self.cpp_info.libs = ["ilda_lib"]
        self.cpp_info.requires = [
            "boost::headers",
            "onetbb::onetbb",
            "opendxa::opendxa",
            "coretoolkit::coretoolkit",
            "structure-identification::structure-identification",
            "common-neighbor-analysis::common-neighbor-analysis",
            "polyhedral-template-matching::polyhedral-template-matching",
            "nlohmann_json::nlohmann_json",
            "spdlog::spdlog",
            "duckdb::duckdb",
        ]
