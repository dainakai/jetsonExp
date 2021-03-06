using CUDA
using CUDA.CUFFT
using StatsBase
using Images
using Spinnaker

function CuGetVector!(vecArray::CuDeviceArray{Float32,3},corArray::CuDeviceArray{Float32,2},gridNum::Int64,corArrSize::Int64,intrSize::Int64)
    gridIdxx = (blockIdx().x-1)*blockDim().x + threadIdx().x
    gridIdxy = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if gridIdxx <= gridNum-1 && gridIdxy <= gridNum-1
        x0::Int64 = 0
        y0::Int64 = 0
        
        tmp::Float32 = 0.0
        @inbounds for j in 1:corArrSize
            for i in 1:corArrSize
                if corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j] > tmp
                    x0 = corArrSize*(gridIdxx-1)+j
                    y0 = corArrSize*(gridIdxy-1)+i
                    tmp = corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j]
                end
            end
        end

        if x0 == 1 || x0 == corArrSize*(gridNum-1) || y0 == 1 || y0 == corArrSize*(gridNum-1)
        # if true
            vecArray[gridIdxy,gridIdxx,1] = Float32(x0) - Float32(intrSize)/2.0 -1.0  - (gridIdxx-1)*corArrSize
            vecArray[gridIdxy,gridIdxx,2] = Float32(y0) - Float32(intrSize)/2.0 -1.0  - (gridIdxy-1)*corArrSize

            return nothing
        else
            valy1x0::Float32 = corArray[y0+1,x0]
            valy0x0::Float32 = corArray[y0,x0]
            valyInv1x0::Float32 = corArray[y0-1,x0]
            valy0x1::Float32 = corArray[y0,x0+1]
            valy0xInv1::Float32 = corArray[y0,x0-1]

            if (valy1x0-2.0*valy0x0+valyInv1x0 == 0.0) || (valy0x1-2.0*valy0x0+valy0xInv1 == 0.0)
                valy0x0 += 0.00001
            end

            vecArray[gridIdxy, gridIdxx,1] = Float32(x0) - (valy0x1 - valy0xInv1)/(valy0x1-2.0*valy0x0+valy0xInv1)/2.0 - Float32(intrSize)/2.0 -1.0  - (gridIdxx-1)*corArrSize
            vecArray[gridIdxy, gridIdxx,2] = Float32(y0) - (valy1x0 - valyInv1x0)/(valy1x0-2.0*valy0x0+valyInv1x0)/2.0 - Float32(intrSize)/2.0 -1.0  - (gridIdxy-1)*corArrSize
            return nothing
        end
    end
end

function CuGetCrossCor!(corArray::CuDeviceArray{Float32,2},img1::CuDeviceArray{Float32,2},img2::CuDeviceArray{Float32,2},gridIdxy::Int64,gridNum::Int64,srchSize::Int64,intrSize::Int64,gridSize::Int64)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if x <= (srchSize-intrSize+1)*(gridNum-1) && y <= (srchSize-intrSize+1)
        gridIdxx = div(x-1,(srchSize-intrSize+1))+1
        idxx = x - (gridIdxx-1)*(srchSize-intrSize+1)
        idxy = y

        a1::Int64 = gridIdxy*gridSize-div(intrSize,2)
        a2::Int64 = gridIdxx*gridSize-div(intrSize,2)
        b1::Int64 = (gridIdxy-1)*gridSize+idxy-1
        b2::Int64 = (gridIdxx-1)*gridSize+idxx-1

        meanA::Float32 = 0.0
        meanB::Float32 = 0.0
        num::Float32 = 0.0
        denomA::Float32 = 0.0
        denomB::Float32 = 0.0

        @inbounds for j in 1:intrSize
            for i in 1:intrSize
                meanA += img1[a1+i,a2+j]
                meanB += img2[b1+i,b2+j]
            end
        end
        meanA /= Float32(intrSize^2)
        meanB /= Float32(intrSize^2)

        @inbounds for j in 1:intrSize
            for i in 1:intrSize
                num += (img1[a1+i,a2+j]-meanA)*(img2[b1+i,b2+j]-meanB)
                denomA += (img1[a1+i,a2+j]-meanA)^2
                denomB += (img2[b1+i,b2+j]-meanB)^2
            end
        end

        corArray[y+(srchSize-intrSize+1)*(gridIdxy-1),x] = num/(CUDA.sqrt(denomA)*CUDA.sqrt(denomB))
    end
    return nothing
end

"""
GPU CUDA Function
???????????? `image1` ??? `image2` ???PIV???????????????????????????????????????1???x, 2???y???
"""
function getPIVMap_GPU(image1::Array{Float32,2}, image2::Array{Float32,2}, imgLen = 1024, gridSize = 128, intrSize = 128, srchSize = 256)
    gridNum = div(imgLen, gridSize)
    corArray = CuArray{Float32}(undef,((srchSize-intrSize+1)*(gridNum-1),(srchSize-intrSize+1)*(gridNum-1)))
    vecArray = CuArray{Float32}(undef,(gridNum-1,gridNum-1,2))

    blockSize = 8
    threads1 = (blockSize,blockSize)
    blocks1 = (cld((srchSize-intrSize+1)*(gridNum-1),blockSize),cld((srchSize-intrSize+1),blockSize))
    blocks2 = (cld(gridNum-1,blockSize),cld(gridNum-1,blockSize))

    d_img1 = cu(image1)
    d_img2 = cu(image2)

    for idx in 1:gridNum-1
        @cuda threads = threads1 blocks = blocks1 CuGetCrossCor!(corArray,d_img1,d_img2,idx,gridNum,srchSize,intrSize,gridSize)
        synchronize()
        # sleep(0.1)
    end

    # display(corArray)

    @cuda threads = threads1 blocks = blocks2 CuGetVector!(vecArray,corArray,gridNum,srchSize-intrSize+1,intrSize)
    synchronize()
    output = Array(vecArray)
    
    return output
end

"""
2??????????????????
"""
function configSetup(camList::CameraList, imgLen = 1024, exposure = 400.0, expratio = 0.5, gain = 0.0)
    cam1 = camList[0]
    cam2 = camList[1]
    
    gammenset1 = gammaenable!(cam1,false)
    @assert gammenset1 == false
    gammenset2 = gammaenable!(cam2,false)
    @assert gammenset2 == false

    adcbits!(cam1, "Bit12")
    @assert adcbits(cam1) == "Bit12" 
    adcbits!(cam2, "Bit12")
    @assert adcbits(cam2) == "Bit12"

    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam1,"AcquisitionFrameRateEnable"),false)
    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam2,"AcquisitionFrameRateEnable"),false)


    pixelformat!(cam1, "Mono16")
    @assert pixelformat(cam1) == "Mono16"
    pixelformat!(cam2, "Mono16")
    @assert pixelformat(cam2) == "Mono16"

    imagedims!(cam1,(imgLen,imgLen))
    imagedims!(cam2,(imgLen,imgLen))
    offsetdims!(cam1,(div(2048-imgLen,2),div(1536-imgLen,2)))
    offsetdims!(cam2,(660,610))
    # offsetdims!(cam2,(div(2048-imgLen,2),div(1536-imgLen,2)))

    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam2,"ReverseY"),true)
    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam2,"ReverseX"),false)

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    exposure!(cam1, exposure)
    exposure!(cam2, exposure*expratio)

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"GainAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"GainAuto"),"Off")
    gain!(cam1, gain)
    gain!(cam2, gain)

    triggermode!(cam1, "On")
    triggersource!(cam1, "Software")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"TriggerSelector"),"FrameStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line1")
    Spinnaker.line_mode(cam1, "Output")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line2")
    Spinnaker.v3_3_enable(cam1,true)
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"Off")

    triggermode!(cam2, "On")
    triggersource!(cam2, "Line3")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerSelector"),"FrameStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"ReadOut")

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"AcquisitionMode"),"Continuous")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"AcquisitionMode"),"Continuous")

    return nothing
end

function loadholo(path::String)
    out::Array{Float32,2} = Float32.(channelview(Gray.(load(path))))
end

function CuTransSqr(datLen::Int64, wavLen::Float32, dx::Float32, Plane::CuDeviceArray{Float32,2})
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = 1.0 - ((x-datLen/2)*wavLen/datLen/dx)^2 - ((y-datLen/2)*wavLen/datLen/dx)^2
    end
    return nothing
end

function CuTransFunc(z0::Float32, wavLen::Float32, datLen::Int64, d_sqrPart::CuDeviceArray{Float32,2}, Plane::CuDeviceArray{ComplexF32,2})
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = exp(2im*pi*(z0)/wavLen*sqrt(d_sqrPart[x,y]))
    end
    return nothing
end

function CuUpdateImposed(datLen::Int64, input::CuDeviceArray{Float32,2},imposed::CuDeviceArray{Float32,2})
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        if input[x,y] < imposed[x,y]
            imposed[x,y] = input[x,y]
        end
    end
    return nothing
end

function fft2d!(plan,out::CuArray{ComplexF32,2}, array::CuArray{Float32,2}, transF::CuArray{ComplexF32,2})
    out .= CUFFT.fftshift(plan*array) .* transF
    return nothing
end

function ifft2d!(plan,array::CuArray{ComplexF32,2},out::CuArray{Float32,2})
    out .= Float32.(abs.(plan*CUFFT.fftshift(array)))
    return nothing
end

function getImposed(img::Array{Float32,2},transF::CuArray{ComplexF32,2},transInt::CuArray{ComplexF32,2},imgLen::Int64=512,blockSize::Int64=16)
    datLen = imgLen*1
    threads = (blockSize,blockSize)
    blocks = (cld(datLen,blockSize),cld(datLen,blockSize))

    # d_img = CuArray{Float32}(undef,(datLen,datLen))
    # plan = CUFFT.plan_fft(d_img)
    # iplan = CUFFT.plan_ifft(d_img)
    # d_img .= 0.5
    # tmpdata = cu(img)
    # d_img .= mean(tmpdata)
    # d_img[div(imgLen,2)+1:div(imgLen,2)+imgLen,div(imgLen,2)+1:div(imgLen,2)+imgLen] .= tmpdata
    d_img = cu(img)
    holo = CuArray{ComplexF32}(undef,(datLen,datLen))
    impImg = CUDA.ones(datLen,datLen)

    # holo .= CUFFT.fftshift(plan*d_img) .* transF
    # fft2d!(plan,holo,d_img,transF)
    holo .= CUFFT.fftshift(CUFFT.fft(d_img)) .* transF
    for idx in 1:100
        holo .= holo .* transInt
        # ifft2d!(iplan,holo,d_img)
        # d_img .= Float32.(abs.(iplan*CUFFT.fftshift(holo)))
        d_img .= Float32.(abs.(CUFFT.ifft(CUFFT.fftshift(holo))))
        @cuda threads = threads blocks = blocks CuUpdateImposed(datLen,d_img,impImg)
    end

    # CUFFT.cufftDestroy(plan)
    # CUFFT.cufftDestroy(iplan)
    output = Array(impImg)
    # return img
    return output
    # return output[div(imgLen,2)+1:div(imgLen,2)+imgLen,div(imgLen,2)+1:div(imgLen,2)+imgLen]
end

# using ImageView
# function testFunc()
#     img1 = loadholo("cam1.bmp")
#     img2 = loadholo("cam2.bmp")
#     imgLen = 1024
#     zF::Float32 = 100.0*1000
#     dz::Float32 = 50.0
#     wavLen::Float32 = 0.532
#     dx::Float32 = 3.45/0.5
#     datLen = imgLen*2
#     blockSize = 16
#     threads = (blockSize,blockSize)
#     blocks = (cld(datLen,blockSize),cld(datLen,blockSize))
#     sqr = CuArray{Float32}(undef,(datLen,datLen))
#     transF = CuArray{ComplexF32}(undef,(datLen,datLen))
#     transInt = CuArray{ComplexF32}(undef,(datLen,datLen))
#     @cuda threads = threads blocks = blocks CuTransSqr(datLen,wavLen,dx,sqr)
#     @cuda threads = threads blocks = blocks CuTransFunc(zF,wavLen,datLen,sqr,transF)
#     @cuda threads = threads blocks = blocks CuTransFunc(dz,wavLen,datLen,sqr,transInt)
#     @time imp1 = getImposed(img1,transF,transInt,imgLen)
#     imshow(imp1)
# end