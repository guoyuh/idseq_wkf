#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang


import os

from DAGflow.dagflow.dag import DAG,Task,ParallelTask

from DAGflow.dagflow import do_dag


inputs = ['1.fasta', "2.fasta", "3.fasta", "4.fasta"]
db = "db.fasta"
db = os.path.abspath(db)

# create a DAG object
my_dag = DAG("blast")
# create the first task 'make_db'
make_db = Task(
    id="make_db",  # your task id, should be unique
    work_dir=".",  # you task work directory
    type="local",  # the way your task run. if "sge", task will submit with qsub
    option={},  # the option of "sge" or "local"
    script="makeblastdb -in %s -dbtype nucl" % db  # the command of the task
)

# when you create a task, then add it to DAG object
my_dag.add_task(make_db)

# then add blast tasks
blast_tasks = ParallelTask(id="blast",
                            work_dir="{id}",
                            type="sge",
                            option="-pe smp 4 -q all.q",
                            script="blastn -in {query} -db %s -outfmt 6 -out {query}.m6",
                            query=inputs)

my_dag.add_task(*blast_tasks)
make_db.set_downstream(*blast_tasks)
# add blast_join task to join blast results
blast_join = Task(
    id="blast_join",
    work_dir=".",
    type="local",  # option is default
    script="cat */*.m6 > blast.all.m6"
)
# you should always remember to add you task to DAG object when created
my_dag.add_task(blast_join)
# this task need a list of tasks in blast_task all done
blast_join.set_upstream(*blast_tasks)

# all of you tasks were added to you workflow, you can run it
do_dag(my_dag)