import cv2
import os

# Get the list of all file names
image_folder = 'Muchang/images2'
images = [img for img in os.listdir(image_folder) if img.endswith(".png")]

# Sort the images by name
images.sort()

# Read the first image to get the width and height
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

# Define the codec and create a VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # or use 'XVID'
video_name = 'output_video2.mp4'
video = cv2.VideoWriter(video_name, fourcc, 50, (width, height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
