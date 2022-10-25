using ClimatePref
using ClimatePref.Units
using Unitful
using Unitful.DefaultSymbols

@test_nowarn create_reference(0.5)
ref = create_reference(0.5)
@test size(ref.array) == (721, 361)
@test ref.array[1:end] == collect(1:260281)