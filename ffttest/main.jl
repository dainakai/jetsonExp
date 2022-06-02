using CUDA,CUDA.CUFFT
using StatsBase
using FFTW
using Random

function fft2d!(array::CuArray{ComplexF32,2},n)
    tmp = CuArray{ComplexF32}(undef,n)
    for i in 1:n
        CUFFT.fft!(array[:,i])
    end

    for i in 1:n
        CUFFT.fft!(array[i,:])
    end
    return nothing
end

function main()
    n = 1024
    x = ComplexF32.(CUDA.rand(n,n))
    @time fftlib = CUFFT.fft(x)
    @time fft2d!(x,n)
    display(mean(abs2.(fftlib.-x)))
end

main()

function testFFTW()
    x = rand(1024,1024)
    x = ComplexF32.(x)
    for i in 1:100
        y = FFTW.fft(x)
    end
end