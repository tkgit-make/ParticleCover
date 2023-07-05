from cv2 import imread, VideoWriter_fourcc, VideoWriter, destroyAllWindows
import os

# Get the list of all file names
image_folder = 'Muchang/images2'
images = [img for img in os.listdir(image_folder) if img.endswith(".png")]

# Sort the images by name
images.sort()

# Read the first image to get the width and height
frame = imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

# Define the codec and create a VideoWriter object
fourcc = VideoWriter_fourcc(*'mp4v') # or use 'XVID'
video_name = 'output_video2.mp4'
video = VideoWriter(video_name, fourcc, 50, (width, height))

for image in images:
    video.write(imread(os.path.join(image_folder, image)))

destroyAllWindows()
video.release()
