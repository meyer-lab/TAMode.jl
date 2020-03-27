using Pkg
Pkg.add("Coverage")
using Coverage

Pkg.activate(".")
Pkg.test("TAMode"; coverage=true)

coverage = process_folder()
LCOV.writefile("coverage-lcov.info", coverage)