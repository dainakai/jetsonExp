using Images, ImageView, Gtk

img = load("./cameraCalib/sample.jpg")
img2 = load("./cameraCalib/sample.jpg")
img3 = load("./cameraCalib/sample.jpg")
# guidict = imshow(img)
# guidict2 = imshow(img3)
# canvas2 = guidict2["gui"]["canvas"]
# canvas = guidict["gui"]["canvas"]

gui = imshow_gui((1500,300),(1,2))
canvases = gui["canvas"]
# anns = [annotations() annotations()]
# annotate!(anns[1],canvases[1,1],AnnotationText(300,0,"Camera 1",fontsize=10))
# idx = annotate!(gui, AnnotationText(300,0,"Camera 1",fontsize=10))

# zr = guidict["roi"]["zoomregion"]
# for j in 1:1024
#     for i in 1:1024
#         img[i,j] = Gray(1.0)
#         # push!(canvas, value(img).currentview
#         imshow(canvas,img)
#     end
# end
# sleep(10)
# img[200:500,200:500] .= 1.0
# imshow(canvas,img)

# while true
#     img[200:500,200:500] .= maximum(img2)
#     sleep(0.00001)
#     imshow(canvas,img)
#     img .= img2
#     sleep(0.00001)
#     imshow(canvas,img)
# end

while true
    img[200:500,200:500] .= maximum(img2)
    img3[200:500,200:500] .= minimum(img2)
    sleep(0.1)
    # imshow(canvas,img)
    # imshow(canvas2,img3)
    gui1 = imshow(canvases[1,1],img)
    # idx = annotate!(gui1, AnnotationText(300,0,"Camera 1",fontsize=10))
    imshow(canvases[1,2],img3)
    Gtk.showall(gui["window"])
    img .= img2
    img3 .= img2
    sleep(0.1)
    # imshow(canvas,img)
    # imshow(canvas2,img3)
    imshow(canvases[1,1],img)
    imshow(canvases[1,2],img3)
    Gtk.showall(gui["window"])
end