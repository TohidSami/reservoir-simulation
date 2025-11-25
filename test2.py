import matplotlib.pyplot as plt
import numpy as np

# Re-creating the Grid class logic internally for demonstration
class Grid:
    def __init__(self, nx=15, ny=15):
        self.NX, self.NY, self.NZ = nx, ny, 1
        self.ACT = [[[1 for _ in range(nx)] for _ in range(ny)] for _ in range(1)]
        self.LA = 0 
        
        mid_y = self.NY // 2
        mid_x = self.NX // 2
        r = 2
        x_start = max(0, mid_x - r)
        x_end = min(self.NX, mid_x + r + 1)
        y_start = max(0, mid_y - r)
        y_end = min(self.NY, mid_y + r + 1)
        
        for j in range(y_start, y_end):
            self.ACT[0][j][mid_x] = -1
        for i in range(x_start, x_end):
            self.ACT[0][mid_y][i] = -1
            
        # Count LA roughly
        for row in self.ACT[0]:
            for val in row:
                if val == 1: self.LA += 1

my_grid = Grid(nx=15, ny=15)
actnum_layer = np.array(my_grid.ACT[0])

plt.figure(figsize=(6, 6))
plt.imshow(actnum_layer, cmap='viridis', origin='lower')
plt.colorbar(label='Activity (1=Active, -1=Inactive)')
plt.title(f'Grid Activity Map ({my_grid.NX}x{my_grid.NY})\nTotal Active Blocks: {my_grid.LA}')
plt.xlabel('X Index')
plt.ylabel('Y Index')
plt.savefig('grid_output.png')

plt.show()