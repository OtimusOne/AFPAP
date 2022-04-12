# File name: AFPAP_p2rank_gallery.py
# Description: Python script for generating P2Rank pocket gallery
# Author: Maghiar Octavian
# Date: 04-04-2022
import argparse
import logging
import pathlib
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
    baseHeight = 480
    hpercent = (baseHeight/float(im.size[1]))
    wsize = int((float(im.size[0])*float(hpercent)))
    im = im.resize((wsize, baseHeight), Image.ANTIALIAS)
    return im


def generatePocketGallery(args):
    df = pd.read_csv(f"{args.outputDir}/work/p2rank_predictions.csv")
    if df.shape[0] == 0:
        return

    pocketImages = ["p2rank_1.png", "p2rank_2.png", "p2rank_3.png",
                    "p2rank_4.png", "p2rank_5.png", "p2rank_6.png",
                    "p2rank_7.png", "p2rank_8.png", "p2rank_9.png", "p2rank_10.png"]
    pocketImagesCrop = ["p2rank_crop_1.png", "p2rank_crop_2.png", "p2rank_crop_3.png",
                        "p2rank_crop_4.png", "p2rank_crop_5.png", "p2rank_crop_6.png", "p2rank_crop_7.png", "p2rank_crop_8.png", "p2rank_crop_9.png", "p2rank_crop_10.png"]

    for img, imgCrop in zip(pocketImages, pocketImagesCrop):
        im = Image.open(f'{args.outputDir}/work/visualizations/{img}')
        im = trim(im)
        im.save(f'{args.outputDir}/work/visualizations/{imgCrop}', "PNG")

    with open(f"{args.AFPAPpath}/config/pocketViewer_template.html", 'r') as template:
        templateData = template.read()
        colors = createSpectrum(df.shape[0], 0.6, 0.6, 1.2)
        pocketDescription = "<span class=\"pocketDescription\">"
        for i in range(df.shape[0]):
            pocketDescription += f"<span class=\"pocketDescriptionSpan\">Pocket{i+1}<span class=\"pocketDescriptionBullet\" style=\"color:rgb({colors[i][0]*255},{colors[i][1]*255},{colors[i][2]*255});\">&#11044</span></span>\t "
        pocketDescription += "</span>"
        templateData = templateData.replace("--pocketColors--", pocketDescription)

        for x, img in enumerate(pocketImagesCrop):
            with open(f'{args.outputDir}/work/visualizations/{img}', "rb") as image_file:
                imgData = base64.b64encode(image_file.read())
                templateData = templateData.replace(f"--pocket{x+1}--", imgData.decode('utf-8'))

        with open(f"{args.outputDir}/work/multiqc_files/p2rank_viewer_mqc.html", 'w') as f:
            f.write(templateData)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")

    args = parser.parse_args()
    consoleLogger = logging.StreamHandler()
    consoleLogger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    fileLogger = logging.FileHandler(f"{args.outputDir}/workflow.log", mode='a')
    fileLogger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)-8.8s - %(message)s",
        handlers=[
            fileLogger,
            consoleLogger
        ]
    )
    logging.info("P2Rank Pocket Viewer...")
    generatePocketGallery(args)


if __name__ == '__main__':
    main()
