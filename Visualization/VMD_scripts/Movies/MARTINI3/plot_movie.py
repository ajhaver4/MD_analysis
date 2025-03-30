import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

# Load CSV Data
data = pd.read_csv(sys.argv[1],comment='#',sep='\t')  # Replace with your CSV file
stepsize = 2
data = data[::stepsize]
x_values = data["Timestep"].values

x_values = np.arange(0, len(x_values), 1)
y_values = data["d1"].values

mask1 = data['RMSD-BtoA'] <= 2.0 
# Define specific indices where scatter points should appear
scatter_indices = np.where(mask1)[0]


scatter_x = x_values[mask1]
scatter_y = y_values[mask1]


# Load video
video_path = sys.argv[2]
cap = cv2.VideoCapture(video_path)

# Set higher resolution
cap.set(cv2.CAP_PROP_FRAME_WIDTH, 1920)  
cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 1080) 

# Get frame rate of the video
fps = int(cap.get(cv2.CAP_PROP_FPS))
frame_interval = 1000 / fps  # Milliseconds per frame

# Get the total number of frames
frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))


print("Number of frames in video: ",frame_count)
frame_count = len(x_values)
print("No. of values in dataframe: ",len(x_values))
# Set up Matplotlib figure with two subplots
fig, ax = plt.subplots(1, 2, figsize=(12, 6), dpi=200)  # Video on left, plot on right
ax_video, ax_plot = ax

# Initialize video frame
ret, frame = cap.read()
if ret:
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)  # Convert BGR to RGB
    img_display = ax_video.imshow(frame)
ax_video.axis("off")
# ax_video.set_title("Video")

# Initialize animated plot
line, = ax_plot.plot([], [], color="navy", lw=2,alpha=0.8)
scatter = ax_plot.scatter([], [], color="black", s=100, marker='*')  # Initially empty scatter plot

ax_plot.set_xlim(min(x_values), max(x_values))
ax_plot.set_ylim(min(y_values), max(y_values))
# ax_plot.set_title("Animated Plot")
fig.suptitle("MARTINI3 - Solution (3D)",fontsize=16)

ax_plot.set_xlabel("Frame", fontsize=15)
ax_plot.set_ylabel(r'$d_1 (nm)$', fontsize=15)
ax_plot.tick_params(axis="both", which="major", labelsize=15)
ax_plot.set_yticks([0,5,10,15,20,25,30])
# Add a frame counter text
# frame_text = ax_plot.text(0.7, 1.04, '', transform=ax_plot.transAxes, fontsize=12, color='black')

fig.tight_layout()
# Animation update function
def update(frame_idx):
    global cap
    ret, frame = cap.read()
    if not ret:
        print("End of video reached. Restarting...")
        print("Frame index: ",frame_idx)
        print(cv2.CAP_PROP_POS_FRAMES)
        cap.set(cv2.CAP_PROP_POS_FRAMES, 0)  # Restart video if it ends
        ret, frame = cap.read()

    # Update video frame
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    img_display.set_array(frame)

    # Update plot
    if frame_idx < len(x_values):
        line.set_data(x_values[:frame_idx], y_values[:frame_idx])

    # Show only scatter points that have been reached
    reached_indices = scatter_indices[scatter_indices <= frame_idx]
    scatter.set_offsets(np.c_[x_values[reached_indices], y_values[reached_indices]])

    # Update frame counter text
    # frame_text.set_text(f'Frame: {frame_idx}/{frame_count}')

    return img_display, line, scatter

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(x_values)+1, interval=frame_interval, blit=False)
# ani = animation.FuncAnimation(fig, update, frames=min(len(x_values), frame_count), interval=frame_interval, blit=False)



# Save animation as an MPEG file
output_video = "Solution_animation.mpeg"
# writer = animation.FFMpegWriter(fps=fps)
# ani.save(output_video, writer=writer)
writer = animation.FFMpegWriter(fps=fps, bitrate=5000)
ani.save("Solution_highres_animation.mp4", writer=writer, dpi=200)
print(f"Animation saved as {output_video}")



# plt.show()

# Release the video capture when done
cap.release()

