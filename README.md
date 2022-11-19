# DEMP-generator-pion-pole-model
Event generator for deep exclusive pion production based on pion pole formula and pion form factor

Event generation of  e p --> e' n pi$^+$  data

# Usage
**Using ROOT to execute the codes:**
>root -l -q test.cpp

**Using ACLiC to compile and execute:**
>root -l

[0] .x test.cpp+

**Using g++ to compile and execute:**
>sed 's/test()/main()/' test.cpp > test2.cpp

>g++ test2.cpp -o run_sim \`root-config --cflags\` \`root-config --libs\`

>./run_sim


# In this version (v1.0), the beam-crossing angle is implemented.
The electron goes to the z direction, consistent with the fixed-target experiment.
see PionExculsiveElectroproduction() or SetBeamCrossAngle(double _angle)
for details of how the crossing angle is implemented.


Download ROOT and get started: https://root.cern.ch/

coding issue contact: rwang@impcas.ac.cn




