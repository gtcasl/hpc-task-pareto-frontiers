#/bin/bash

echo "threads,time,average_power,max_power" > fair.csv
for x in `seq 1 1 112`
do 
    f=cholesky.fair.$x.log
    t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> fair.csv
done

echo "threads,time,average_power,max_power" > coreconstrained.csv
for x in `seq 2 1 112`
do 
    f=cholesky.coreconstrained.$x.log
    t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> coreconstrained.csv
done

echo "limit,time,average_power,max_power" > poweraware.csv
for x in `seq 30 5 300`
do 
    f=cholesky.poweraware.$x.log
    t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> poweraware.csv
done

echo "limit,time,average_power,max_power" > slackaware.csv
for x in `seq 30 5 300`
do 
    f=cholesky.slackaware.$x.log
    t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m" >> slackaware.csv
done
