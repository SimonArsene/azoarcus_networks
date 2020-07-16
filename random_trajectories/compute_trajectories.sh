nb_trajectories=1000
nb_points=20
list_spec="4 8 12 16 20 24"
for nb_spec in $list_spec 
do 
touch network_trajectories_alph_${nb_spec}.txt 
done
date
python compute_trajectories.py $nb_trajectories $nb_points -L $list_spec
date