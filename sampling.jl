import TAMode
using Test
using Turing

f = open("log.txt", "w")
redirect_stdout(f)
redirect_stderr(f)

@testset "Test sampling." begin
    samp = sample(TAMode.A549model, NUTS(0.65), 500)
    write("chain-file-10_1.jls", samp)
end

close(f)
