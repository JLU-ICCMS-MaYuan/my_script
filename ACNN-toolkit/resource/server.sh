#!/bin/bash

# remote server setting
# Ensure that server SSH public key authentication passes
set_server() {
    case "$1" in
        server1)    # for cpu server
            remote_ip="hostname@remote_ip"
            remote_port="22"
            task_sys="slurm"     # slurm, pbs
            remote_work_place="/remote/workplace/label_plat"
            priv_rsa_path="~/.ssh/id_rsa"
            ;;
        server2)    # gpu server for training
            remote_ip="hostname@remote_ip"
            remote_port="22"
            task_sys="slurm"     # slurm, pbs
            remote_work_place="/remote/workplace/label_plat"
            priv_rsa_path="~/.ssh/id_rsa"
            ;;
        *)
            echo "Invalid server"
            ;;
    esac
}

unset_server() {
    unset remote_ip remote_port remote_work_place task_sys priv_rsa_path
}

echo_server() {
    echo_line1
    echo "Current Server Info:"
    echo "Remote IP: $remote_ip"
    echo "Remote Port: $remote_port"
    echo "Task System: $task_sys"
    echo "Remote Work Place: $remote_work_place"
    echo "Private RSA Path: $priv_rsa_path"
    echo_line1
}
