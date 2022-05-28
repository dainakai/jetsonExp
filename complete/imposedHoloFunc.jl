using CUDA

function CuTransSqr(datLen, wavLen, dx, Plane)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = 1.0 - ((x-datLen/2)*wavLen/datLen/dx)^2 - ((y-datLen/2)*wavLen/datLen/dx)^2
    end
    return nothing
end

function CuTransFunc(z0, wavLen, datLen, d_sqrPart, Plane)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = exp(2im*pi*(z0)/wavLen*sqrt(d_sqrPart[x,y]))
    end
    return nothing
end

function CuUpdateImposed(datLen, input,imposed)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        if input[x,y] < imposed[x,y]
            imposed[x,y] = input[x,y]
        end
    end
    return nothing
end

function getImposed(img,datLen=1024,wavLen=0.6328,dx=10.0,zF=120.0*1000,dz=50.0)
    blockSize = 16
    threads = (blockSize,blockSize)
    blocks = (cld(datLen,blockSize),cld(datLen,blockSize))

    d_img = cu(img)
    sqr = CuArray{Float32}(undef,(datLen,datLen))
    transF = CuArray{ComplexF32}(undef,(datLen,datLen))
    transInt = CuArray{ComplexF32}(undef,(datLen,datLen))
    holo = CuArray{ComplexF32}(undef,(datLen,datLen))
    impImg = CUDA.ones(datLen,datLen)

    @cuda threads = threads blocks = blocks CuTransSqr(datLen,wavLen,dx,sqr)
    @cuda threads = threads blocks = blocks CuTransFunc(zF,wavLen,datLen,sqr,transF)
    @cuda threads = threads blocks = blocks CuTransFunc(dz,wavLen,datLen,sqr,transInt)

    holo = CUFFT.fftshift(CUFFT.fft(d_img)) .* transF
    for idx in 1:1000
        holo = holo .* transInt
        d_img = Float32.(abs.(ifft(fftshift(holo))))
        @cuda threads = threads blocks = blocks CuUpdateImposed(datLen,d_img,impImg)
        # save(string("./reconstInv/",lpad(idx,5,"0"),".bmp"),img)
    end

    return Array(impImg)
end