#!/bin/bash 

function SIGINT_handler {                                                          
    echo -e "\nSIGINT Caught ..."                                                        
    #Propagate the signal up to the shell                                       
    kill -s SIGINT $$                                                           
    # 130 is the exit status from Ctrl-C/SIGINT                                 
    exit 130                                                                    
}                                                                               
                                                                                
function exit_handler {                                                         
    echo "... Script Exiting!"                                                       
}                                                                               
                                                                                
checkExitCode() {                                                               
  exitvalue=$?
  if [ ! $exitvalue -eq 0 ]; then                                                       
    error_exit " The last function return a wrong exit code, execution halted. Check logs for further details."
  fi                                                                            
}                                                                               
                                                                                
error_exit() {                                                                  
  echo -e "$(date): ${1:-" Unknown Error, execution halted "}" 2>&1 | tee -a doitgromacs_err.log 
  exit 1                                                                        
}                              

trap 'SIGINT_handler; exit 0' INT                                                  
#trap 'checkExitCode; exit 1' ERR                                                                                                 
trap 'exit_handler' EXIT     

