using CairoMakie
using Spinnaker
using Images
using Gtk, Gtk.ShortNames

"""
2カメラの設定
"""
function configSetup(camList, exposure = 200.0, gain = 0.0)
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

    # Spinnaker.set!(Spinnaker.SpinBooleanNode(cam2,"ReverseX"),true)

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    exposure!(cam1, exposure)
    exposure!(cam2, exposure)

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

function CuGetVector!(vecArray::CuDeviceArray{Float32,3},corArray::CuDeviceArray{Float32,2},gridNum::Int64,corArrSize::Int64,intrSize::Int64)
    gridIdxx = (blockIdx().x-1)*blockDim().x + threadIdx().x
    gridIdxy = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if gridIdxx <= gridNum-1 && gridIdxy <= gridNum-1
        x0::Int64 = 0
        y0::Int64 = 0
        
        tmp::Float32 = 0.0
        for i in 1:corArrSize
            for j in 1:corArrSize
                if corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j] > tmp
                    x0 = corArrSize*(gridIdxx-1)+j
                    y0 = corArrSize*(gridIdxy-1)+i
                    tmp = corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j]
                end
            end
        end

        valy1x0::Float32 = corArray[y0+1,x0]
        valy0x0::Float32 = corArray[y0,x0]
        valyInv1x0::Float32 = corArray[y0-1,x0]
        valy0x1::Float32 = corArray[y0,x0+1]
        valy0xInv1::Float32 = corArray[y0,x0-1]

        if (valy1x0-2.0*valy0x0+valyInv1x0 == 0.0) || (valy0x1-2.0*valy0x0+valy0xInv1 == 0.0)
            valy0x0 += 0.00001
        end

        vecArray[gridIdxy, gridIdxx,1] = Float32(x0) - (valy0x1 - valy0xInv1)/(valy0x1-2.0*valy0x0+valy0xInv1)/2.0 - Float32(intrSize)/2.0 -1.0  - (gridIdxx-1)*corArrSize
        vecArray[gridIdxy, gridIdxx,2] = Float32(y0) - (valy1x0 - valyInv1x0)/(valy1x0-2.0*valy0x0+valyInv1x0)/2.0 - intrSize/2.0 -1.0  - (gridIdxy-1)*corArrSize
    end
    return nothing
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

        for i in 1:intrSize
            for j in 1:intrSize
                meanA += img1[a1+i,a2+j]
                meanB += img2[b1+i,b2+j]
            end
        end
        meanA /= Float32(intrSize^2)
        meanB /= Float32(intrSize^2)

        for i in 1:intrSize
            for j in 1:intrSize
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
画像配列 `image1` と `image2` のPIVマップを計算し、返します。1はx, 2はy。
"""
function getPIVMap_GPU(image1, image2, imgLen = 1024, gridSize = 128, intrSize = 128, srchSize = 256)
    gridNum = div(imgLen, gridSize)
    corArray = CuArray{Float32}(undef,((srchSize-intrSize+1)*(gridNum-1),(srchSize-intrSize+1)*(gridNum-1)))
    vecArray = CuArray{Float32}(undef,(gridNum-1,gridNum-1,2))

    blockSize = 16
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

function imgAcquisition(camList,state,stopbutton)
    cam1 = camList[0]
    cam2 = camList[1]

    start!(cam1)
    start!(cam2)
    trigger!(cam1)
    img1 = getimage(cam1)
    img2 = getimage(cam2)
    arr1 = CameraImage(img1,Float32, normalize = true)
    arr2 = CameraImage(img2,Float32, normalize = true)
    vecArray = getPIVMap_GPU(arr1,arr2)

    f = Figure(resolution = (1600,500),figure_padding = 1)
    arrowax = Makie.Axis(f[1,3], aspect = 1 , yreversed=false, backgroundcolor="white")
    imageax1 = Makie.Axis(f[1, 1], aspect = DataAspect(), yreversed = false, title = "Camera 1")
    imageax2 = Makie.Axis(f[1, 2], aspect = DataAspect(), yreversed = false, title = "Camera 2")
    imgObservable1 = Observable(rotr90(RGB.(arr1,arr1,arr1)))
    imgObservable2 = Observable(rotr90(RGB.(arr2,arr2,arr2)))
    image!(imageax1,imgObservable1)
    image!(imageax2,imgObservable2)

    vecxObservable = Observable(rotr90(vecArray[:,:,1]))
    vecyObservable = Observable(-rotr90(vecArray[:,:,2]))
    strObservable = Observable(vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2)))

    xs = [i*128 for i in 1:n]
    ys = [i*128 for i in 1:n]
    arrows!(arrowax, xs,ys, vecxObservable,vecyObservable, arrowsize=10, lengthscale=20, arrowcolor = strObservable, linecolor = strObservable)

    Makie.save("./loopfig.pdf",f)
    
    signal_connect(x -> state = "stop", stopbutton, "clicked")
    while state == "running"
        # while state == "running"
        trigger!(cam1)
        getimage!(cam1,img1)
        getimage!(cam2,img2)
        arr1 = CameraImage(img1,Float32, normalize = true)
        arr2 = CameraImage(img2,Float32, normalize = true)
        vecArray = getPIVMap_GPU(arr1,arr2)
        imgObservable1[] = rotr90(RGB.(arr1,arr1,arr1))
        imgObservable2[] = rotr90(RGB.(arr2,arr2,arr2))
        vecxObservable[] = rotr90(vecArray[:,:,1])
        vecyObservable[] = -rotr90(vecArray[:,:,2])
        strObservable[] = vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2))
        sleep(0.01)
        # i += 1
        # sleep(0.01)
    end
    stop!(cam1)
    stop!(cam2)
    Spinnaker._release!(cam1)
    Spinnaker._release!(cam2)
    # close(f)
end

function main()
    # get camera list
    if !(@isdefined camList)
        camList = CameraList()
        display(camList)
        configSetup(camList)
    end
    configSetup(camList)

    win = Window("Camera Controller")
    v = GtkBox(:v)
    l1 = GtkLabel("")
    b1 = GtkButton("Start")
    l2 = GtkLabel("")
    b2 = GtkButton("Stop")
    l3 = GtkLabel("")
    b3 = GtkButton("Finish")
    l4 = GtkLabel("")
    push!(win, v)
    push!(v,l1)
    push!(v, b1)
    push!(v,l2)
    push!(v, b2)
    push!(v,l3)
    push!(v, b3)
    push!(v,l4)

    showall(win)

    imgAcquisTask(camList, state, stopbutton) = @task begin; imgAcquisition(camList,state,stopbutton); end

    function start(camList, stopbutton)
        state = "running"
        schedule(imgAcquisTask(camList, state, stopbutton))
        println("Start clicked!")
    end

    signal_connect(x -> start(camList, b2), b1, "clicked")
    signal_connect(x -> exit(), b3, "clicked")

    if !isinteractive()
        c = Condition()
        signal_connect(win, :destroy) do widget
            notify(c)
        end
        wait(c)
    end
end

main()