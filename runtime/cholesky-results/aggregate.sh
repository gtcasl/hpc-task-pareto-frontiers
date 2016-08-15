#/bin/bash

echo "threads,time,average_power,max_power,estimated_max_power" > fair.csv
for x in `seq 1 1 112`
do 
    f=cholesky.fair.$x.log
    t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    e=`tail -n 5 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
    echo "$x,$t,$a,$m,$e" >> fair.csv
done

echo "threads,time,average_power,max_power,estimated_max_power" > coreconstrained.csv
for x in `seq 1 1 112`
do 
    f=cholesky.coreconstrained.$x.log
    if [ -f $f ]
    then
        t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        e=`tail -n 5 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        echo "$x,$t,$a,$m,$e" >> coreconstrained.csv
    fi
done

echo "limit,time,average_power,max_power,estimated_max_power" > poweraware.csv
for x in `seq 30 5 330`
do 
    f=cholesky.poweraware.$x.log
    if [ -f $f ]
    then
        t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        e=`tail -n 5 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        echo "$x,$t,$a,$m,$e" >> poweraware.csv
    fi
done

echo "limit,time,average_power,max_power,estimated_max_power" > slackaware.csv
for x in `seq 30 5 330`
do 
    f=cholesky.slackaware.$x.log
    if [ -f $f ]
    then
        t=`tail -n 2 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        a=`tail -n 4 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        m=`tail -n 3 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        e=`tail -n 5 $f | head -n 1 | awk -F "[=,]" '{print $2}'`
        echo "$x,$t,$a,$m,$e" >> slackaware.csv
    fi
done
