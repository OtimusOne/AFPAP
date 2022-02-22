import base64
from PIL import Image, ImageChops


def trim(im):
    bg = Image.new(im.mode, im.size, im.getpixel((0, 0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        im = im.crop(bbox)
    baseHeight = 360
    hpercent = (baseHeight/float(im.size[1]))
    wsize = int((float(im.size[0])*float(hpercent)))
    im = im.resize((wsize, baseHeight), Image.ANTIALIAS)
    return im


imgPath = ["./output/work/p2rank_1.png", "./output/work/p2rank_2.png", "./output/work/p2rank_3.png",
           "./output/work/p2rank_4.png", "./output/work/p2rank_5.png", "./output/work/p2rank_6.png"]
imgCropPath = ["./output/work/p2rank_crop_1.png", "./output/work/p2rank_crop_2.png", "./output/work/p2rank_crop_3.png",
               "./output/work/p2rank_crop_4.png", "./output/work/p2rank_crop_5.png", "./output/work/p2rank_crop_6.png"]

for path, cropPath in zip(imgPath, imgCropPath):
    im = Image.open(path)
    im = trim(im)
    im.save(cropPath, "PNG")


with open('config/pocketViewer_template.html', 'r') as template:
    templateData = template.read()
    for x, img in enumerate(imgCropPath):
        with open(img, "rb") as image_file:
            imgData = base64.b64encode(image_file.read())
            templateData = templateData.replace(
                f"--pocket{x+1}--", imgData.decode('utf-8'))

    with open('./output/work/p2rank_viewer_mqc.html', 'w') as f:
        f.write(templateData)
