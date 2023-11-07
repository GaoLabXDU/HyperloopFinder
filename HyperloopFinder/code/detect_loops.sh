chrom=$1
loop_dir=$2
resolution=$3
assembly=$4
downweighting=$5
cluster_file=$6
max_cluster_size=$7
loop_lowerbound=$8
loop_upperbound=${9}
fithic_pass=${10}
filter_by_loop=${11}
correct_bias=${12}

echo $chrom
echo $loop_dir
echo $resolution


contacts=${loop_dir}/${chrom}_contacts_${resolution}.tsv.gz
fragments=${loop_dir}/${chrom}_fragments_${resolution}.tsv.gz
bias=${loop_dir}/${chrom}_bias_${resolution}.tsv.gz


python  get_sprite_contacts.py   --assembly	${assembly}\
                                 --chromosome	${chrom}\
                                 --resolution	${resolution}\
                                 --downweighting	${downweighting}\
                                 --clusters	${cluster_file}\
                                 --max_cluster_size	${max_cluster_size}\
                                 --raw_contacts	${contacts}\
                                 --fragments ${fragments}

if [ $filter_by_loop -eq 1 ] || [ $correct_bias -eq 1 ]; then
       python ../utils/HiCKRy.py  -i ${contacts}\
                            -f ${fragments}\
                            -o ${bias}\
                            -x 0.05

fi

if [ $filter_by_loop -eq 1 ]; then
       fithic -i ${contacts}\
              -f ${fragments}\
              -t ${bias}\
              -o ${loop_dir}/${chrom}_fithic\
              -r ${resolution}\
              -p $fithic_pass\
              -v \
              -x intraOnly\
              -L $loop_lowerbound\
              -U $loop_upperbound
fi
