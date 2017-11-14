# Lightweaver 

<p>
Light weaver is a ray tracer built to test out some code for 3dda acceleration <br/>
3dda is a 3d modification of dda (digital differential analyzer) which is used to quickly determine which pixels on a 2d plane are intersected by a line. This in 3 dimensions can be used to quickly traverse an octtree reducing the amount of triangle intersection tests needed. 
<p/>

## Build

mkdir build <br/> 
cd build <br/>
cmake ../src <br/>
make -j8 <br/>

## Run 

./LightWeaver <filename> <br/>
./LightWeaver ../data/monkey.txt


## Example renders

![alt text](https://i.imgur.com/pHzup87.png "Stairs")

![alt text](https://i.imgur.com/D26wuiq.png "Starship")

![alt text](https://i.imgur.com/1L8Mert.png "Monkey")

![alt text](https://i.imgur.com/EObguxH.png "Xwing")


