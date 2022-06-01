using CairoMakie
using Spinnaker
using Images
using Gtk, Gtk.ShortNames
using CUDA

include("imposedHoloFunc.jl")

function imgAcquisition(camList,state,stopbutton,imgLen,gridSize,intrSize, srchSize)
    # Gabor Reconstruction Prepareing
    zF::Float32 = 100.0*1000
    dz::Float32 = 50.0
    wavLen::Float32 = 0.532
    dx::Float32 = 3.45/0.5
    datLen = imgLen*2
    blockSize = 16
    threads = (blockSize,blockSize)
    blocks = (cld(datLen,blockSize),cld(datLen,blockSize))
    sqr = CuArray{Float32}(undef,(datLen,datLen))
    transF = CuArray{ComplexF32}(undef,(datLen,datLen))
    transInt = CuArray{ComplexF32}(undef,(datLen,datLen))
    @cuda threads = threads blocks = blocks CuTransSqr(datLen,wavLen,dx,sqr)
    @cuda threads = threads blocks = blocks CuTransFunc(zF,wavLen,datLen,sqr,transF)
    @cuda threads = threads blocks = blocks CuTransFunc(dz,wavLen,datLen,sqr,transInt)
    println("Gabor Init OK")
    ############ end ###############


    cam1 = camList[0]
    cam2 = camList[1]

    start!(cam1)
    start!(cam2)
    trigger!(cam1)
    img1 = getimage(cam1)
    img2 = getimage(cam2)
    stop!(cam1)
    stop!(cam2)
    arr1 = Array(CameraImage(img1,Float32, normalize = true))
    arr2 = Array(CameraImage(img2,Float32, normalize = true))
    display(size(arr1))
    vecArray = getPIVMap_GPU(arr1,arr2,imgLen,gridSize,intrSize,srchSize)
    println("PIV OK")

    f = Figure(resolution = (1600,500),figure_padding = 1)
    imageax1 = Makie.Axis(f[1, 1], aspect = DataAspect(), yreversed = false, title = "Camera 1")
    imageax2 = Makie.Axis(f[1, 2], aspect = DataAspect(), yreversed = false, title = "Camera 2")
    arrowax = Makie.Axis(f[1,3], aspect = 1 , yreversed=false, backgroundcolor="white")
    imgObservable1 = Observable(rotr90(RGB.(arr1,arr1,arr1)))
    imgObservable2 = Observable(rotr90(RGB.(arr2,arr2,arr2)))
    image!(imageax1,imgObservable1)
    image!(imageax2,imgObservable2)

    vecxObservable = Observable(rotr90(@view vecArray[:,:,1]))
    vecyObservable = Observable(-rotr90(@view vecArray[:,:,2]))
    strObservable = Observable(vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2)))

    n = div(imgLen,gridSize) -1
    xs = [i*gridSize for i in 1:n]
    ys = [i*gridSize for i in 1:n]
    arrows!(arrowax, xs,ys, vecxObservable,vecyObservable, arrowsize=10, lengthscale=20, arrowcolor = strObservable, linecolor = strObservable)
    Makie.limits!(arrowax,0,gridSize*(n+1),0,gridSize*(n+1))

    println("Makie Output OK")
    Makie.save("./loopfig.pdf",f)
    
    signal_connect(x -> state = "stop", stopbutton, "clicked")
    while state == "running"
        # while state == "running"
        display(rand())
        start!(cam1)
        start!(cam2)
        trigger!(cam1)
        getimage!(cam1,img1)
        getimage!(cam2,img2)
        stop!(cam1)
        stop!(cam2)
        arr1 = Array(CameraImage(img1,Float32, normalize = true))
        arr2 = Array(CameraImage(img2,Float32, normalize = true))
        # imp1 = arr1
        imp1 = getImposed(arr1,transF,transInt,imgLen,blockSize)
        # imp2 = arr2
        imp2 = getImposed(arr2,transF,transInt,imgLen,blockSize)
        vecArray = getPIVMap_GPU(imp1,imp2,imgLen,gridSize,intrSize,srchSize)
        imgObservable1[] = rotr90(RGB.(imp1,imp1,imp1))
        imgObservable2[] = rotr90(RGB.(imp2,imp2,imp2))
        vecxObservable[] = rotr90(vecArray[:,:,1])
        vecyObservable[] = -rotr90(vecArray[:,:,2])
        strObservable[] = vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2))
        sleep(0.01)
        Makie.save("./loopfig.pdf",f)
    end
    Spinnaker._release!(cam1)
    Spinnaker._release!(cam2)
    println("Stop Clicked!")
end

function main()
    imglen = 512
    gridsize = div(imglen,8)
    intrsize = div(imglen,8)
    srchsize = div(imglen,4)

    # get camera list
    if !(@isdefined camList)
        camList = CameraList()
        display(camList)
        configSetup(camList,imglen)
    end

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

    imgAcquisTask(camList, state, stopbutton) = @task begin; imgAcquisition(camList,state,stopbutton,imglen,gridsize,intrsize,srchsize); end

    function start(camList, stopbutton)
        state = "running"
        schedule(imgAcquisTask(camList, state, stopbutton))
        println("Start clicked!")
    end

    signal_connect(x -> start(camList, b2), b1, "clicked")
    signal_connect(x -> exit(), b3, "clicked")

    showall(win)

    if !isinteractive()
        c = Condition()
        signal_connect(win, :destroy) do widget
            notify(c)
        end
        wait(c)
    end
end

main()