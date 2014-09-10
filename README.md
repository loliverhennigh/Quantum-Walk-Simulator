# Quantum Walk Simulator

This program is designed to make a visual quantum walk and allow the user to test graphs for perfect or pretty good state transfer. When the program runs it will ask the user to input a file name. This file contains the graph the user wants to run the walk on. You are allowed real and imaginary weights on the graph (our research was on imaginary weighted graphs). Note that the adjacency matrix for the graph must be hermitian (to keep it unitary). Here is an example of a file (with comments).



## Two Node Graph

Here is what a graph file lookes like

```
2 (number of nodes)

(Adjacency Matrix)

-1  0  1.414  0  (first two values are the real and then imag part of a number) 

1.414  0  -1  0  

(This matrix is 2 by 2)

(position of node in space) 

0 0 0

1 0 0

(This is used my the graphics side)
 
(starting vector)

1 0 0 0  (again, first two numbers for each value, real, imag)
```

When running this file it produces this screen

![ScreenShot](https://github.com/loliverhennigh/Quantum-Walk-Simulator/blob/master/two_node.png)

Here is another example of a file called pres.txt

![ScreenShot](https://github.com/loliverhennigh/Quantum-Walk-Simulator/blob/master/pres.png)


## How it works

This algorithm uses spectral decomposition for its speed and accuracy. Taylor approximation is to slow and diverges quickly. Here is a brief look at the math.


```
|distribution> = e^(-itA) dot |start>
```

Where, 

```
A is the adjacency matrix, 

t is time

|start> (vector) is the initial distribution. In the example above its (1,0) so it starts on the first node

|distribution> is the final distribution on the nodes. In the example above at time roughly 1.11 we see perfect state transfer to the second 
not, ie (0,1).
```


## Commands

Here is a list of commands

```
t = timestep of .001

f = timestep of .01

g = timestep of .1

- = reverse direction

1 = fullscreen

2 = exit fullscreen

arrows = rotate

w,a,s,d = move the image
```


## How to compile




Here are the steps,

You will need freeglut and gsl. To install these on ubunto simply type

```
sudo apt-get install freeglut3-dev
```

and

```
sudo apt-get install libgsl0-dev
```

Next compile the program

```
make
```

and to exacute

```
./walks
```

If glut gives you any trouble consult http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/









