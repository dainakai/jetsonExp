using Images, ImageView, Gtk, Gtk.ShortNames
using Spinnaker

```
SpinEnumNode SpinFloatNode SpinBooleanNode SpinIntegerNode
CameraNodeMap CameraTLStreamNodeMap CameraTLDeviceNodeMap
```

# get camera list
camList = CameraList()
println(camList)

function configSetup(camList, exposure = 6001.0, gain = 13.10)
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
    # Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"TriggerSelector"),"AcquisitionStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line1")
    Spinnaker.line_mode(cam1, "Output")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line2")
    Spinnaker.v3_3_enable(cam1,true)
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"Off")

    triggermode!(cam2, "On")
    triggersource!(cam2, "Line3")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerSelector"),"FrameStart")
    # Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerSelector"),"AcquisitionStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"ReadOut")

    # Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"AcquisitionMode"),"SingleFrame")
    # Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"AcquisitionMode"),"SingleFrame")

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"AcquisitionMode"),"Continuous")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"AcquisitionMode"),"Continuous")
end

function imgAcquisition(camList,state, button)
    gui = imshow_gui((1500,300),(1,2))
    canvases = gui["canvas"]

    cam1 = camList[0]
    cam2 = camList[1]

    # reset!(cam1)
    # reset!(cam2)

    # try
    #     while state == "running"
    #         start!(cam1)
    #         start!(cam2)
    #         trigger!(cam1)
    #         img1 = getimage(cam1)
    #         img2 = getimage(cam2)
    #         arr1 = CameraImage(img1,Float64, normalize = true)
    #         arr2 = CameraImage(img2,Float64, normalize = true)
    #         println("img1 maximum: ",maximum(arr1))
    #         println("img2 maximum: ",maximum(arr2))
    #         println("")
    #         arr1[1,1] = 0.99
    #         arr2[1,1] = 0.99
    #         # println("Camera 1 timestamp: ",timestamp(arr1))
    #         # println("Camera 2 timestamp: ",timestamp(arr2))
    #         imshow(canvases[1,1], arr1)
    #         imshow(canvases[1,2], arr2)
    #         sleep(0.01)
    #         Gtk.showall(gui["window"])
    #         stop!(cam1)
    #         stop!(cam2)
    #     end
    # catch
    #     stop!(cam1)
    #     stop!(cam2)
    #     ImageView.closeall()
    # end
    signal_connect(x -> state = "stop", button, "clicked")
    while state == "running"
        start!(cam1)
        start!(cam2)
        trigger!(cam1)
        img1 = getimage(cam1)
        img2 = getimage(cam2)
        arr1 = CameraImage(img1,Float64, normalize = true)
        arr2 = CameraImage(img2,Float64, normalize = true)
        println("img1 maximum: ",maximum(arr1))
        println("img2 maximum: ",maximum(arr2))
        println("")
        arr1[1,1] = 0.99
        arr2[1,1] = 0.99
        # println("Camera 1 timestamp: ",timestamp(arr1))
        # println("Camera 2 timestamp: ",timestamp(arr2))
        imshow(canvases[1,1], arr1)
        imshow(canvases[1,2], arr2)
        sleep(0.01)
        Gtk.showall(gui["window"])
        stop!(cam1)
        stop!(cam2)
        println("stop ok")
    end
    # ImageView.closeall()
    # reset!(cam1)
    # reset!(cam2)
end

win = Window("Count Click")
v = Box(:v)
l = Label("")
# b = Button("Start")
# l2 = Label("")
b = Button("Stop")
l3 = Label("")
push!(win, v)
push!(v,l)
# push!(v, b)
# push!(v,l2)
push!(v, b)
push!(v,l3)
# set_gtk_property!(v, :expand, b, true)

showall(win)

# signal_connect(x -> imgAcquisition(camList,state), b, "clicked")
# signal_connect(x -> state = "stop", b2, "clicked")

configSetup(camList)

# function start()
    
# end

# imgAcquisition(camList,"running")
state = "running"
imgAcquisition(camList, state, b)

if !isinteractive()
    c = Condition()
    signal_connect(win, :destroy) do widget
        notify(c)
    end
    wait(c)
end