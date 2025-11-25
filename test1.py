import matplotlib.pyplot as plt
import numpy as np
from mygrid import Grid

# Instantiate the Grid
my_grid = Grid(nx=15, ny=15)

# Extract ACT (Activity) for plotting (Layer 0)

actnum_layer = np.array(my_grid.ACT[0])

# Plotting
plt.figure(figsize=(6, 6))
plt.imshow(actnum_layer, cmap='jet', origin='lower') 
plt.colorbar(label='Activity (1=Active, -1=Inactive)')
plt.title(f'Grid Activity Map ({my_grid.NX}x{my_grid.NY})\nTotal Active Blocks (LA): {my_grid.LA}')
plt.xlabel('X Index')
plt.ylabel('Y Index')
plt.grid(which='major', color='white', linestyle='-', linewidth=0.5)
plt.xticks(np.arange(-0.5, 15, 1)) 
plt.yticks(np.arange(-0.5, 15, 1)) 
# Save the figure
plt.savefig('grid_activity_check.png')
print("Plot saved as grid_activity_check.png")

plt.show()