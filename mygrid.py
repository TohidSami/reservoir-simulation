class Grid:
    def __init__(self,nx,ny, nz):
        self.NX,self.NY,self.NZ=nx,ny ,nz
        self.dx , self.dy, self.dz=30, 30 ,10
        #m to ft 
        self.dx=3.28084*self.dx
        self.dy=3.28084*self.dy
        self.dz=3.28084*self.dz

        self.POR=0.1
        self.PERM_Y_all=100
        self.PERM_X_all=100
        self.PERM_Z_all=100

        self.NT = self.NY * self.NX * self.NZ
        self.LA=0
        self.Vb = self.dx * self.dy * self.dz
        self.C_DARCY = 0.00112712
        self.ACT=self.matrix3D()
        # 2. Initialize Arrays (Empty Structures)
        self.PROS = self.matrix3D()
        self.PERM_X = self.matrix3D()
        self.PERM_Y = self.matrix3D()
        self.PERM_Z = self.matrix3D()
        
        self.IT_IJ = self.matrix3D()
        self.IA_IJ = self.matrix3D()

        self.IJ_IT = []
        self.IJ_IA = []
        self.NL = []
        self.TL = []
        self.TRANS = []
        self.PV = []

        self.properties()
        self.connectivity()


    def matrix3D(self):
            return [[[1 for _ in range(self.NX)] for _ in range(self.NY)] for _ in range(self.NZ)]
    def matrix2D(self,row,col):
            return [[1 for j in range(col)] for i in range(row)]

    def properties(self):
        for z in range(self.NZ):
                for j in range(self.NY):
                    for i in range(self.NX):
                        self.PROS[z][j][i] = self.POR
                        self.PERM_X[z][j][i] = self.PERM_X_all
                        self.PERM_Y[z][j][i] = self.PERM_Y_all
                        self.PERM_Z[z][j][i] = self.PERM_Z_all
                        self.ACT[z][j][i] = 1
        mid_y=self.NY//2
        mid_x=self.NX//2
        r=2
        x_start=max(0 , mid_x-r)
        x_end=min(self.NX,mid_x+r)
        y_start=max(0 , mid_y-r)
        y_end=min(self.NY,mid_y+r)
        for j in range(y_start,y_end+1):
            for z in range(self.NZ):
                self.ACT[z][j][mid_x]=-1
        for i in range(x_start,x_end+1):
            for z in range(self.NZ):
                self.ACT[z][mid_y][i]=-1
        
        
    def connectivity(self):    
        self.IJ_IT = self.matrix2D(self.NT, 3)
        IT = 0
        IA = 0
        # temp_IA_IJ = self.matrix3D()
        for z in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    IT +=1
                    self.IT_IJ[z][j][i]=IT
                    self.IJ_IT[IT-1][0]=i
                    self.IJ_IT[IT-1][1]=j
                    self.IJ_IT[IT-1][2]=z
                    if self.ACT[z][j][i]>=0:
                        IA +=1
                        self.IA_IJ[z][j][i] = IA
                        # temp_IA_IJ[z][j][i] = IA
                    else:
                        self.IA_IJ[z][j][i] = -1
                        # temp_IA_IJ[z][j][i] = IA                        
        self.LA=IA
        self.IJ_IA=self.matrix2D(self.LA,3)
        self.NL=self.matrix2D(self.LA,6)
        self.TL=self.matrix2D(self.LA,6)
        self.PV=[[0.0] for _ in range(self.LA)]
        self.TRANS=self.matrix2D(self.LA,6)
        ac=0
        for z in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    if self.ACT[z][j][i] > 0:
                         self.IJ_IA[ac][0]=i
                         self.IJ_IA[ac][1]=j
                         self.IJ_IA[ac][2]=z
                         self.PV[ac][0]=self.NX*self.NY*self.NZ*self.POR
        
        
                         if (i-1 < 0) or (self.ACT[z][j][i-1]<0):self.NL[ac][0]=-1
                         else: self.NL[ac][0]=self.IA_IJ[z][j][i-1] 
                        
                         if (i+1 > self.NX-1) or (self.ACT[z][j][i+1] < 0): self.NL[ac][1] = -1
                         else: self.NL[ac][1] = self.IA_IJ[z][j][i+1]
                        
                         if (j-1 < 0) or (self.ACT[z][j-1][i] < 0): self.NL[ac][2] = -1
                         else: self.NL[ac][2] = self.IA_IJ[z][j-1][i]
                                               
                         if (j+1 > self.NY-1) or (self.ACT[z][j-1][i] < 0): self.NL[ac][3] = -1
                         else: self.NL[ac][3] = self.IA_IJ[z][j+1][i]

                         if (z-1 < 0) or (self.ACT[z-1][j][i] < 0): self.NL[ac][4] = -1
                         else: self.NL[ac][4] = self.IA_IJ[z-1][j][i]

                         if (z+1 > self.NZ-1) or (self.ACT[z+1][j][i] < 0): self.NL[ac][5] = -1
                         else: self.NL[ac][5] = self.IA_IJ[z+1][j][i]
                         for k in range(6):
                          if self.NL[ac][k] != -1: 
                              self.TL[ac][k]=self.NL[ac][k]
                          else: self.TL[ac][k]=-1
                         ac +=1
        for o in range(IA):
            ix=self.IJ_IA[o][0]
            iy=self.IJ_IA[o][1]
            iz=self.IJ_IA[o][2]
            kx=self.PERM_X[iz][iy][ix]
            ky=self.PERM_Y[iz][iy][ix]
            kz=self.PERM_Z[iz][iy][ix]
            tx=(self.C_DARCY*kx*self.dy*self.dz)/(self.dx/2)
            ty=(self.C_DARCY*ky*self.dx*self.dz)/(self.dy/2)
            tz=(self.C_DARCY*kz*self.dx*self.dy)/(self.dz/2)
            for k in range(6):
                if k==0 or k==1:
                    if self.NL[o][k]>0:
                        self.TL[o][k]=tx
                if k==2 or k==3:
                    if self.NL[o][k]>0:
                        self.TL[o][k]=ty
                if k==4 or k==5:
                    if self.NL[o][k]>0:
                        self.TL[o][k]=tz
        for o in range(IA):
            for k in range(6):
                # o = o+1
                if k==0:
                  right=self.NL[o][0]-1
                  if right>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[right][1]
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))
                  else: self.TRANS[o][k]=0
                if k==1:
                  left=self.NL[o][k]-1
                  if left>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[left][0]
                      if t2!=t1:
                          print("slam")
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))      
                  else: self.TRANS[o][k]=0
                if k==2:
                  front=self.NL[o][k]-1
                  if front>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[front][3]
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))
                  else: self.TRANS[o][k]=0
                if k==3:
                  back=self.NL[o][k]-1
                  if back>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[back][2]
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))
                  else: self.TRANS[o][k]=0
                if k==4:
                  up=self.NL[o][k]-1
                  if up>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[up][5]
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))                        
                  else: self.TRANS[o][k]=0
                if k==5:
                  down=self.NL[o][k]-1
                  if down>=0:
                      t1=self.TL[o][k]
                      t2=self.TL[down][4]
                      self.TRANS[o][k] = 1 / ((1/t1) + (1/t2))
                  else: self.TRANS[o][k]=0        
        print("dkufh")     


