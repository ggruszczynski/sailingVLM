# pyLLTandVLM


[![CircleCI](https://circleci.com/gh/ggruszczynski/sailingVLM/tree/main.svg?style=svg)](https://circleci.com/gh/ggruszczynski/sailingVLM/?branch=main)

python implementation of a 3D Vortex Lattice Method

```
            The geometry is described using the following CSYS.
            lem - leading edge (luff) of the main sail
            tem - trailing edge (leech) of the main sail

                        Z ^ (mast)     
                          |
                         /|
                        / |
                       /  |              ^ Y     
                   lem_NW +--+tem_NE    / 
                     /    |   \        /
                    /     |    \      /
                   /      |     \    /
                  /       |      \  /
                 /        |       \/
                /         |       /\
               /          |      /  \
              /           |     /    \
             /            |    /      \
            /      lem_SW |---/--------+tem_SE
           /              |  /          
  (bow) ------------------|-/-------------------------| (stern)
         \                |/                          |
    ------\---------------*---------------------------|-------------------------> X (water level)
        
```
