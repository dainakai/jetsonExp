using GLMakie
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

function imgAcquisition(camList,state,stopbutton)
    cam1 = camList[0]
    cam2 = camList[1]

    start!(cam1)
    start!(cam2)
    trigger!(cam1)
    img1 = getimage(cam1)
    # println("img1 ok")
    img2 = getimage(cam2)
    # println("img2 ok")
    arr1 = CameraImage(img1,Float32, normalize = true)
    arr2 = CameraImage(img2,Float32, normalize = true)
    arr1 = RGB.(arr1,arr1,arr1)
    arr2 = RGB.(arr2,arr2,arr2)
    # display(arr1)
    # display(arr2)
    # arr1[1,1] = 0.99
    # arr2[1,1] = 0.99
    # arr1[1,2] = 0.0
    # arr2[1,2] = 0.0
    # stop!(cam1)
    # stop!(cam2)
    # Spinnaker._release!(cam1)
    # Spinnaker._release!(cam2)

    f = Figure(resolution = (1000,500),figure_padding = 1)

    imageax1 = Makie.Axis(f[1, 1], aspect = DataAspect(), yreversed = false, title = "Camera 1")
    imageax2 = Makie.Axis(f[1, 2], aspect = DataAspect(), yreversed = false, title = "Camera 2")
    imgObservable1 = Observable(rotr90(arr1))
    imgObservable2 = Observable(rotr90(arr2))
    image!(imageax1,imgObservable1)
    image!(imageax2,imgObservable2)

    display(f)
    # i = 0
    signal_connect(x -> state = "stop", stopbutton, "clicked")
    while state == "running"
        # while state == "running"
        trigger!(cam1)
        getimage!(cam1,img1)
        getimage!(cam2,img2)
        arr1 = CameraImage(img1,Float32, normalize = true)
        arr2 = CameraImage(img2,Float32, normalize = true)
        arr1 = RGB.(arr1,arr1,arr1)
        println(rand())
        arr2 = RGB.(arr2,arr2,arr2)
        imgObservable1[] = rotr90(arr1)
        imgObservable2[] = rotr90(arr2)
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