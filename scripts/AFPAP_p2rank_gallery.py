import base64
from turtle import color
from PIL import Image, ImageChops
import colorsys
import pandas as pd


def createSpectrum(numColors,  sat,  hueMin,  hueMax):
    table = [0]*numColors

    if (numColors == 1):

        table[0] = createSpectrum(2, sat, hueMin, hueMax)[0]

    else:
        hueRange = hueMax - hueMin

        for i in range(numColors):
            hue = hueMin + hueRange * (i / numColors)

            if (hue > 1):
                hue -= 1

            if (i % 2 == 0):
                table[i] = colorsys.hsv_to_rgb(hue, sat,  0.9)
            else:
                table[i] = colorsys.hsv_to_rgb(hue, sat,  0.7)
    return table


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


def generatePocketGallery():

    df = pd.read_csv("./output/work/p2rank_predictions.csv")

    if df.shape[0] == 0:
        return

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
        colors = createSpectrum(df.shape[0], 0.6, 0.6, 1.2)
        pocketDescription = "<span style=\"width: 100%;display: inline-flex;justify-content: center;\">"
        for i in range(df.shape[0]):
            pocketDescription += f"<span style=\"padding: 0 5px;\">Pocket{i+1}<span style=\"color:rgb({colors[i][0]*255},{colors[i][1]*255},{colors[i][2]*255});font-size:16px\">&#11044</span></span>\t "
        pocketDescription += "</span>"
        templateData = templateData.replace(
            "--pocketColors--", pocketDescription)

        for x, img in enumerate(imgCropPath):
            with open(img, "rb") as image_file:
                imgData = base64.b64encode(image_file.read())
                templateData = templateData.replace(
                    f"--pocket{x+1}--", imgData.decode('utf-8'))

        with open('./output/work/p2rank_viewer_mqc.html', 'w') as f:
            f.write(templateData)


generatePocketGallery()
