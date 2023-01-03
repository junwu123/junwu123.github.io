Enclosed are two models for the paper [Wu2016]. (The interior of) the kitten model is optimized for stiffness, and (the interior of) the letter P model is optimized for static stability of standing. Please refer to Fig. 7 in the paper.

- In each group, you will find an input mesh which is essentially the outer surface. (-> KittenInput_100.0%.obj)
- From this mesh, an adaptive rhombic structure is initialized. (-> KittenInitial_54.5%.obj)
- Starting from the initial structure, the rhombic cells are further adaptive refined based on an optimization of stiffness or static stability. (-> KittenOptimized_70.7%.obj)
- To reveal the interior structure, the carved rhombic cells are displayed separately. (-> KittenOptimized_Removed_(100-70.7)%.obj)
- The percentage value indicates the volume percentage with respect to the volume enclosed by the outer surface.

[Wu2016] Jun Wu, Charlie C.L. Wang, Xiaoting Zhang, Rüdiger Westermann, Self-supporting rhombic infill structures for additive manufacturing, Computer-Aided Design, Volume 80, November 2016, Pages 32-42, ISSN 0010-4485, http://dx.doi.org/10.1016/j.cad.2016.07.006.

The group of kitten models is based on the kitten model provided courtesy of the AIM@SHAPE Shape Repository.

Oct. 20, 2016, Jun Wu (email: j.wu-1@tudelft.nl)
