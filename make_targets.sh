# generate science targets

python3 target_generator.py -o data/test_sci.dat -r 100000 -s 10000 -a 0. -d 0. --select-single-cobra


# generate sky targets

python3 target_generator.py --plot -o data/test_sky.dat -r 1000 -s 1000 -a 0. -d 0. --select-single-cobra


# generate calibration targets

python3 target_generator.py --plot -o data/test_cal.dat -r 100 -s 1 -a 0. -d 0. --select-single-cobra


# run netflow

python3 netflow_asiaa.py


