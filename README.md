cProj
=====
This project is about image segmentation using parallel graph partitioning algorithm for gray scale
images.
Image segmentation refers to the decomposition of a scene into its components. Segmentation
subdivides an image into its constituent regions or objects. Image segmentation is used in many
applications such as locate objects in satellite images (roads, forests, etc.), face recognition, fingerprint
recognition, traffic control systems, brake light detection, machine vision, etc.
There are many methods for image segmentation such as watershed,Mean-shift,Edge detection &
Graph-Based partitioning. Graph-Based partitioning treats image as undirected weighted graph which
consists of
1. vertexes which are the pixels in feature space.
2.
Edges between each pair of vertexes
3. Weight for each edge represents the similarity between pixels.
The algorithm used in this project is based on graph partitioning methods. The used technique is
minimum spanning tree method but in parallelization way.
Parallelization helps most in making the segmentation algorithm more efficient by reducing the
calculation time which can affects the performance in sequential algorithms especially with large sized
images.
