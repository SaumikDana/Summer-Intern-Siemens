from dolfin import*
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import math
import numpy as np
import sys
import xlrd
from timeit import default_timer as timer

#read in input file
start=timer() #timer start
xval=[]
yval=[]
time=[]
power=[]
workbook=xlrd.open_workbook('data1.xlsx',on_demand=True)
worksheet=workbook.sheet_by_index(0)

count=20

#count=1
#while worksheet.cell_value(count,3)==worksheet.cell_value(count+1,3):
#    count=count+1

for i in range(count): #count is number of data points on the first microslice
    xval.append(worksheet.cell_value(i+1,1)) 
    yval.append(worksheet.cell_value(i+1,2))
    time.append(worksheet.cell_value(i+1,0))
    power.append(worksheet.cell_value(i+1,8))

#get minimum values for further processing
xvalmin=min(xval)
yvalmin=min(yval)
timevalmin=min(time)

for i in range(count):
    xval[i]=xval[i]-xvalmin
    yval[i]=yval[i]-yvalmin
    #assuming data for time given in milliseconds
    time[i]=(time[i]-timevalmin)*10**-3

end=timer() #timer end
print "Preprocessing time=%15.8f"%(end-start)

#homogenized problem
nx=50
ny=50

domain_x=[-0.5*max(xval),1.5*max(xval)]
domain_y=[-0.5*max(yval),1.5*max(yval)]

hx_homo=(domain_x[1]-domain_x[0])/nx
hy_homo=(domain_y[1]-domain_y[0])/ny

mesh=RectangleMesh(Point(domain_x[0],domain_y[0]),Point(domain_x[1],domain_y[1]),nx,ny,'crossed')
V=FunctionSpace(mesh,"CG",1)

#original problem
nx_orig=200
ny_orig=200

hx=(domain_x[1]-domain_x[0])/nx_orig
hy=(domain_y[1]-domain_y[0])/ny_orig

mesh_orig=RectangleMesh(Point(domain_x[0],domain_y[0]),Point(domain_x[1],domain_y[1]),nx_orig,ny_orig,'crossed')

#contour plots for temperature
triang=tri.Triangulation(*mesh.coordinates().reshape((-1,2)).T,triangles=mesh.cells())

u=TrialFunction(V) #trial space for the homogenized problem
v=TestFunction(V) #test space for the homogenized problem

#initial condition
u_i=interpolate(Expression("300",degree=1),V)

#convection coefficient in W/mm^2 K
hc=0.00001

#dirichlet BC
bc=DirichletBC(V,Expression("300",degree=1),DomainBoundary())

sys.setrecursionlimit(20000) #hack for memory issues

#time stepping loop
T=max(time)
t=0
counter=0
u_sol=u_i #initialize finite element solution for the homogenized problem at t=0

k_exp=Expression("0",degree=0) #conductivity units W/mm-K
c_exp=Expression("0",degree=0) #specific heat units J/mm^3-K

#toolpath
xplot=[]
yplot=[]
xplot.append(xval[0])
yplot.append(yval[0])

while counter<count-1: #count-1 is the number of segments

    if(power[counter]==0): #no computation required if laser does nothing
        #segment bounds
        xbound=xval[counter+1]-xval[counter]
        ybound=yval[counter+1]-yval[counter]
        seglength=sqrt(xbound**2+ybound**2)
        t+=dt #increment time
        counter+=1 #increment segment counter
        end=timer() #end time
        #print out everything relevant
        print "Time(sec),Wall clock time elapsed,segment length=%15.8f:%15.8f:%15.8f"%(t,end-start,seglength)
        continue

    dt=time[counter+1]-time[counter] #adaptive time step 

    start=timer() #start time

    xv1=min(xval[counter],xval[counter+1])
    xv2=max(xval[counter],xval[counter+1])
    yv1=min(yval[counter],yval[counter+1])
    yv2=max(yval[counter],yval[counter+1])

    f_exp=Expression('(x[0]>=xv1 & x[0]<=xv2 & x[1]>=yv1-hy & x[1]<=yv2+hy)?f:g',xv1=xv1,xv2=xv2,yv1=yv1,yv2=yv2,hx=hx,hy=hy,f=power[counter],g=0,degree=0)
    #source term

    k_exp_add=Expression('(x[0]>=xv1 & x[0]<=xv2 & x[1]>=yv1-hy & x[1]<=yv2+hy)?f:g',xv1=xv1,xv2=xv2,yv1=yv1,yv2=yv2,hx=hx,hy=hy,f=0.0215,g=0.0003,degree=0)
    #conductivity
    k_exp=Expression('k_exp+k_exp_add',k_exp=k_exp,k_exp_add=k_exp_add,degree=2)

    c_exp_add=Expression('(x[0]>=xv1 & x[0]<=xv2 & x[1]>=yv1-hy & x[1]<=yv2+hy)?f:g',xv1=xv1,xv2=xv2,yv1=yv1,yv2=yv2,hx=hx,hy=hy,f=0.00595,g=0.00425,degree=0)
    #specific heat
    c_exp=Expression('c_exp+c_exp_add',c_exp=c_exp,c_exp_add=c_exp_add,degree=2)

    #homogenized conductivity
    k_eff_1=assemble((1/k_exp)*dx(mesh_orig)) #arrive at integral of conductivity inverse
    k_eff=assemble((1/k_eff_1)*dx(mesh_orig)) #arrive at integral of the inverse of the above

    #RHS for homogenized problem
    L=c_exp*u_sol*v*dx(mesh)+f_exp*dt*v*dx(mesh)
    #contribution from boundary condition
    L=L+dt*hc*u_i*v*ds(mesh)

    #LHS for homogenized problem    
    a=c_exp*u*v*dx(mesh)+k_eff*dt*inner(grad(u),grad(v))*dx(mesh)
    #contribution from boundary condition
    a=a+hc*u*v*ds(mesh)

    #solve homogenized problem
    u_sol=Function(V)

    problem=LinearVariationalProblem(a,L,u_sol)
    solver=LinearVariationalSolver(problem)
    #gmres with ilu preconditioner
    solver.parameters["linear_solver"]="gmres"
    solver.parameters["preconditioner"]="ilu"
    gmres_prm=solver.parameters["krylov_solver"]
    gmres_prm["absolute_tolerance"]=1E-8
    gmres_prm["relative_tolerance"]=1E-4
    gmres_prm["maximum_iterations"]=1000
    solver.solve()

    #append toolpath
    xplot.append(xval[counter+1])
    yplot.append(yval[counter+1])

    #segment bounds
    xbound=xval[counter+1]-xval[counter]
    ybound=yval[counter+1]-yval[counter]
    seglength=sqrt(xbound**2+ybound**2)

    t+=dt #increment time
    counter+=1 #increment segment counter
    end=timer() #end time

    #print out everything relevant
    print "Time(sec),Wall clock time elapsed,segment length,power=%15.8f:%15.8f:%15.8f:%15.8f"%(t,end-start,seglength,power[counter-1])

    #contour plot of temperature
    Z=u_sol.compute_vertex_values(mesh)
    plt.figure
    plt.tricontourf(triang,Z)
    plt.colorbar()

    #overlay plot of toolpath
    plt.plot(xplot,yplot,'-ro',markersize=2)
    plt.show()


