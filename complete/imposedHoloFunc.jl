using CUDA
using CUDA.CUFFT
using StatsBase

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

# function meanPadding(out,in,datLen,mean)
#     x = (blockIdx().x-1)*blockDim().x + threadIdx().x
#     y = (blockIdx().y-1)*blockDim().y + threadIdx().y
#     if x <= div(datLen,4)*3 && y <= div(datLen,4)*3 && x >= div(datLen,2)+1 && y >= div(datLen,2)+1
#         if input[x,y] < imposed[x,y]
#             imposed[x,y] = input[x,y]
#         end
#     end
#     return nothing
# end

function getImposed(img,imgLen=1024,wavLen=0.532,dx=3.45/0.5,zF=100.0*1000,dz=50.0)
    datLen = imgLen*2
    blockSize = 16
    threads = (blockSize,blockSize)
    blocks = (cld(datLen,blockSize),cld(datLen,blockSize))

    tmpdata = cu(img)
    d_img = CuArray{Float32}(undef,(datLen,datLen))
    d_img .= mean(tmpdata)
    d_img[div(imgLen,2)+1:div(imgLen,2)+imgLen,div(imgLen,2)+1:div(imgLen,2)+imgLen] .= tmpdata[:,:]
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
        d_img = Float32.(abs.(CUFFT.ifft(CUFFT.fftshift(holo))))
        @cuda threads = threads blocks = blocks CuUpdateImposed(datLen,d_img,impImg)
        # save(string("./reconstInv/",lpad(idx,5,"0"),".bmp"),img)
    end

    output = Array(impImg)
    return output[div(imgLen,2)+1:div(imgLen,2)+imgLen,div(imgLen,2)+1:div(imgLen,2)+imgLen]
end