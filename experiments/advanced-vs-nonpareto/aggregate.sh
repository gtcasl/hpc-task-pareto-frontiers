#/bin/bash

echo "threads,time,average_power,max_power" > control.all
for x in `seq 1 1 112`
do 
    t=`grep time_s cholesky.control.$x.log | awk -F "[=,]" '{print $2}'`
    a=`grep average_power cholesky.control.$x.log | awk -F "[=,]" '{print $2}'`
    m=`grep max_power cholesky.control.$x.log | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> control.all
done

echo "limit,time,average_power,max_power" > experimental.all
for x in `seq 30 5 220`
do 
    t=`grep time_s cholesky.experimental.$x.log | awk -F "[=,]" '{print $2}'`
    a=`grep average_power cholesky.experimental.$x.log | awk -F "[=,]" '{print $2}'`
    m=`grep max_power cholesky.experimental.$x.log | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> experimental.all
done
