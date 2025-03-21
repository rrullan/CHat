#!/bin/bash
transition_file=transition.txt
grid_size=80
save_cube='No'
cutoff=0.00


echo "                                                       .-.        .-."
echo "                                                     -  )      .'  /"
echo "                                          __        /  /      /   /"
echo "                .--.                   .-'  ''''-._/__(_____.'   /"
echo "             .-'  _)          __     .'  ___.--\'                J"
echo "            /  '-.,   .------\'  '--.J---\'                      /"
echo "            F      '-<                                          .\'"
echo "            '.        '.                                       .--\'"
echo "            .\'|        \\\\                                   . /'"
echo "          .\'            J                                      /"
echo "         /              /                                    _/"
echo "        /               \'                           -\'  -"
echo "      .\'                                          ("
echo "     /                               __.----\'      \\"
echo "    /                         __.---\'       '.      L"
echo " .-\'     _.              .--\'\'                \\     \\"
echo "/       \'             .-\'                      \\     \\"
echo "'--.       /'  __.---\'                          \\     \\"
echo "    \\     /  .\'                                  \\     \\"
echo "    -'.    _/                                     \\     \\"
echo "     \' \\_,-\'-._                                  \\     L"
echo "           \\'.'-                                    '.    L"
echo '                                                      \\   |'
echo '                                                       L   L'
echo '                                                       |   )'
echo "                                                       \\  /"
echo "                                                         '' "


echo "      Welcome to the CHaracterization and Analysis of the electronic Transitions software         "
echo ""
echo "                                    ____ _   _       _   "
echo "                                   / ___| |_| | __ _| |_ "
echo "                                  | |   | '_  |/ _  | __|"
echo "                                  | |___| | | | (_| | |_ "
echo "                                   \____|_| |_|\__'_|\__|"
echo ""
echo " This software has been written by Antonin Dufour (M2 Student) and Raphael Rullan (PhD student)   "
echo " under the supervision of Pr. Tangui le Bahers and Dr. Stephan Steinmann of the LCH ENS de Lyon.  "
echo " The objective of this software is to compute different indexes (Tozer, DCT ...) to analyse the   "
echo " electronic transitions either in periodic systems (CP2K software) or in non-periodic (Orca) ones."
echo ""
echo "The python script uses the Cube-Toolz package: https://github.com/funkymunkycool/Cube-Toolz for"
echo "cube manipulations."
echo "The implementation of the Tozer index is based on the paper J. Chem. Phys. 128, 044118 (2008)."
echo "The implementation of the DCT index is based on the paper J. Chem. Theory Comput. 2011, 7, 8, 2498–2506"
echo ""
echo 'Miaou !'


if [ -f ${transition_file} ] ; then
        echo -n "" > ${transition_file}
fi
index_label="A B a b c 1 2 3 4 5 11 12 13 14 15 16 17 20 21"
function not_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    LIST_WHITESPACES=`echo $LIST | tr "$DELIMITER" " "`
    for x in $LIST_WHITESPACES; do
        if [ "$x" = "$VALUE" ]; then
            return 1
        fi
    done
    return 0
}

while [[ ! $Index_Type =~ 21 ]]
do
        if [ -f ${transition_file} ] ; then
                echo -n "" > ${transition_file}
        fi
        
        echo "    =============================================================="
        echo "         Choose a Code:"
        echo ""
        echo "               - CP2K: A"
        echo "               - Orca: B"
        echo ""
        echo "               - Exit: 21"
        echo "    --------------------------------------------------------------"
        read -p "    Enter a type: " Index_Type
       
        if not_in_list "$index_label"  " "  $Index_Type ; then
                echo ""
                echo "#########################################################################"
                echo ""
                echo "                 ERROR : Please choose a valid option "
                echo ""
                echo "#########################################################################"
                echo ""
               
	fi

	if [ $Index_Type = A ] ; then
		echo ""
		echo "    =============================================================="
		echo ""
		echo "                  Welcome to the CP2K interface"
		echo ""
		echo "    =============================================================="
		echo ""
		if [ $# -eq 0 ] ; then
			
       			echo "    Please enter the name of the INP-file of the CP2K calculation"
        		read -p "    Enter the filename: " aux_file
        		echo "    Please enter the name of the OUT-file of the CP2K calculation"
        		read -p "    Enter the filename: " input_file
        		echo "    Please enter the name of the orbital (cube) obtained with the CP2K calculation"
        		read -p "    Enter the filename: " orb_file
		fi
		while [[ ! ($Index_Type -eq 20) && ! ($Index_Type -eq 21)  ]]
		do
       			echo "    ============================================================="
        		echo "         Grid Size: ${grid_size}x${grid_size}x${grid_size}"
        		echo "         Name of the INP-file: ${aux_file}"
			echo "         Name of the OUT-file: ${input_file}"
			echo "         Name of the orb-file: ${orb_file}"
        		echo "         Cutoff in the transition decomposition: ${cutoff}"
        		echo "    =============================================================="
        		echo "         Type of Index :"
        		echo ""
        		echo "               - Tozer Index: 1"
        		echo "               - DCT Index: 2"
        		echo ""
        		echo ""
			echo "               - Change OUT-file: 11"
			echo "               - Change INP-file: 12"
			echo "               - Change ORB-file: 13"	
			echo "               - Change cutoff: 14"
			echo "               - Name of the Tozer file: 15"
        		echo ""
                        echo "               - Change the code: 20"
                        echo "               - Exit: 21"
        		echo "    --------------------------------------------------------------"
        		read -p "    Enter a type: " Index_Type
        		echo "    =============================================================="
        		if not_in_list "$index_label"  " "  $Index_Type ; then
                		echo ""
                		echo "#########################################################################"
                		echo ""
                		echo "           ERROR : Please choose a valid option "
                		echo ""
                		echo "#########################################################################"
                		echo ""
                	
				fi

				if [ $Index_Type = 11 ] ; then
                			while [[ ! $whiling =~ ^(n|y)$ ]]
                			do
                        			echo "    The default name of the OUT-file would be : ${input_file}"
                        			read -p "    Is this the one you want (y/n) : " whiling
                			done
                			if [ $whiling = n ] ; then
                        			read -p "    Enter the filename : " input_file
                			fi
        			fi

        			if [ $Index_Type = 12 ] ; then
                			whiling=''
                			while [[ ! $whiling =~ ^(n|y)$ ]]
                			do
                        			echo "    The default name of the INP-file would be : ${aux_file}"
                        			read -p "    Is this the one you want (y/n) : " whiling
                			done
                			if [ $whiling = n ] ; then
                        			read -p "    Enter the filename : " aux_file
                			fi
        			fi

				if [ $Index_Type = 13 ] ; then
                                        whiling=''
                                        while [[ ! $whiling =~ ^(n|y)$ ]]
                                        do
                                                echo "    The default name of the ORB-file would be : ${aux_file}"
                                                read -p "    Is this the one you want (y/n) : " whiling
                                        done
                                        if [ $whiling = n ] ; then
                                                read -p "    Enter the filename : " orb_file
                                        fi
                                fi 

        			if [ $Index_Type = 14 ] ; then
                			echo "    The cutoff for the transition decomposition is : ${cutoff}"
                			read -p "    Please enter the new cutoff : " cutoff
				fi

				if [ $Index_Type = 15 ] ; then
                                	echo "    The name for the Tozer file is by default : Tozer-state"
                                        read -p "    Please enter the new suffixe : " name
                                fi

					
				if [ $Index_Type = 1 ]; then
                			read -p "    Excited State Number: " State_Index
                			echo "    =============================================================="
                			echo ""
					sed -i 's/nstates/NSTATES/Ig' ${input_file}
                			nb_state=(`python -c"import excited_index; excited_index.find_nb_state_calc_CP2K('${aux_file}')"`)
                			max_index=${nb_state[0]}
                			

                			if [ ${State_Index} -gt $max_index ] ; then
                        			echo " ------- ERROR ---------"
                        			echo ""
                        			echo "There is no State ${State_Index} in the file ${input_file}"
                        			echo ""
                        			echo " ----------------------"
                        			echo ""
                        		
                			fi

                			if [ ${State_Index} -gt ${nb_state[0]} ] ; then
                        			echo ""	
                			else
                        			sed -i "s/(alp)//g" ${input_file}
                        			list_orb_tmp=$(grep 'WAVEFUNCTION' *.cube > orb.txt)
                        			list_orb=$(awk '{print $3}' orb.txt)
                        			num=$((${State_Index}+1))
                        			get_state=$(grep -A ${num} "number   energy (eV)       x           y           z     strength (a.u.)" ${input_file}|tail -n 1 > tmp.txt)
                        			state_index=$(awk '{print $2}' tmp.txt)
                        			energy_state=$(awk '{print $3}' tmp.txt)
                        			i=1
                        			temp=$(grep -A $i "${state_index}   ${energy_state} eV" ${input_file} | tail -n 1)
                        			echo $temp > tmp.txt
                        			index_orb=$(awk '{print $1}' tmp.txt)
                        			check=$(awk '{print $3}' tmp.txt)
                        			while [[ "$check" != "eV" ]];
                        			do
                                			echo $temp >> ${transition_file}
                                			orb_occ=$(awk '{print $1}' tmp.txt)
                                			coeff_occ=$(awk '{print $3}' tmp.txt)
                               				abs_coeff_occ=$(echo "sqrt($coeff_occ*$coeff_occ)" | bc -l)
                              				cond=$(echo "$abs_coeff_occ > $cutoff" | bc -l)

                                			if  not_in_list "$list_orb"  " "  $orb_occ  ; then
                                       				if [ $cond -gt 0 ] ; then
                                        				echo ""
                                      					echo "#########################################################################"
                                        				echo ""
                                      					echo "           ERROR: The orbital $orb_occ is not found. "
                                       					echo ""
                                        				echo "#########################################################################"
                                        				echo ""
                                        				
                                				fi
                                			fi

                                			orb_vac=$(awk '{print $2}' tmp.txt)

                                			if not_in_list "$list_orb"  " "  $orb_vac ; then
                                        			if [ $cond -gt 0 ] ; then
                                       					echo ""
                                        				echo "#########################################################################"
                                        				echo ""
                                        				echo "           ERROR: The orbital $orb_vac is not found. "
                                        				echo ""
                                        				echo "#########################################################################"
                                        				echo ""
                                         			
                                				fi
                                			fi
                                			((i++))
                                			temp=$(grep -A $i "${state_index}   ${energy_state} eV" ${input_file}| tail -n 1)
                                			echo $temp > tmp.txt
                                			index_orb=$(awk '{print $1}' tmp.txt)
                                			check=$(awk '{print $3}' tmp.txt)

                        			done
                			fi

					python excited_index.py tozer ${orb_file} ${cutoff} CP2K
               				rm  square*
                			rm  mult*.cube
					rm transition.txt
					rm orb.txt
					rm tmp.txt
                			mv Tozer.txt Tozer_state_CP2K_${State_Index}.txt
					rm -r __pycache__

					if [ -n "$name" ]; then
						mv Tozer_state_CP2K_${State_Index}.txt Tozer_state_CP2K_${State_Index}_${name}.txt
					fi

                			elif [ $Index_Type = 2 ]; then
              					echo "     Not implemented yet.    "
                			fi
		done


	
		
	elif [ $Index_Type = B ] ; then
		echo ""
		echo "    =============================================================="
                echo ""
                echo "                  Welcome to the Orca interface"
                echo ""
                echo "    =============================================================="
		echo ""
		if [ $# -eq 0 ] ; then
        	
        		echo "    Please enter the name of the OUT-file of the ORCA calculation"
        		read -p "    Enter the filename: " input_file
		else
        		input_file=${1}
		fi
		filename="${input_file%.*}"
		gbw_file="${filename}.gbw"
		cis_file="${filename}.cis"
		while [[ ! ($Index_Type -eq 20) && ! ($Index_Type -eq 21)  ]]
		do
        		if [ -f ${transition_file} ] ; then
                		echo -n "" > ${transition_file}
        		fi
        		echo ""
        		echo "    ============================================================="
        		echo "         GBW-file name: ${gbw_file}"
        		echo "         CIS-file name: ${cis_file}"
			echo "         OUT-file name: ${input_file}"
        		echo "         Grid Size: ${grid_size}x${grid_size}x${grid_size}"
        		echo "         Auto-Save of CUBE-file: ${save_cube}"
        		echo "         Cutoff in the transition decomposition: ${cutoff}"
        		echo "    =============================================================="
        		echo "         Type of Index:"
        		echo ""
        		echo "               - Tozer Index: 1"
        		echo "               - DCT Index: 2"
        		echo "               - Omega Index For ICS: 3"
        		echo "               - Overlapping of Transitions: 4"
        		echo ""
			echo "               - Only generate the orbitals: 5"
			echo ""
        		echo "               - Change GBW-file: 11"
        		echo "               - Change CIS-file: 12"
			echo "               - Change OUT-file: 13"
        		echo "               - Change Grid Size: 14"
        		echo "               - Auto-Save of CUBE-file: 15"
        		echo "               - Change the cutoff: 16"
			echo "               - Name of the Tozer file: 17"
        		echo ""
        		echo "               - Change the code: 20"
			echo "               - Exit: 21"
        		echo "    --------------------------------------------------------------"
        		read -p "    Enter a type: " Index_Type
        		echo "    =============================================================="
        		if not_in_list "$index_label"  " "  $Index_Type ; then
                		echo ""
                		echo "#########################################################################"
                		echo ""
                		echo "                 ERROR: Please choose a valid option "
                		echo ""
                		echo "#########################################################################"
                		echo ""	
        		fi

        		if [ $Index_Type = 11 ] ; then
                		while [[ ! $whiling =~ ^(n|y)$ ]]
                		do
                        		echo "    The default name of the GBW-file would be : ${gbw_file}"
                        		read -p "    Is this the one you want (y/n) : " whiling
                		done
                		if [ $whiling = n ] ; then
                        		read -p "    Enter the filename : " gbw_file
                		fi
        		fi

        		if [ $Index_Type = 12 ] ; then
                		whiling=''
                		while [[ ! $whiling =~ ^(n|y)$ ]]
                		do
                        		echo "    The default name of the CIS-file would be : ${cis_file}"
                        		read -p "    Is this the one you want (y/n) : " whiling
                		done
                		if [ $whiling = n ] ; then
                        		read -p "    Enter the filename : " cis_file
                		fi
        		fi

			if [ $Index_Type = 13 ] ; then
                                whiling=''
                                while [[ ! $whiling =~ ^(n|y)$ ]]
                                do
                                        echo "    The default name of the OUT-file would be : ${input_file}"
                                        read -p "    Is this the one you want (y/n) : " whiling
                                done
                                if [ $whiling = n ] ; then
                                        read -p "    Enter the filename : " input_file
					filename="${input_file%.*}"
                			gbw_file="${filename}.gbw"
                			cis_file="${filename}.cis"
                                fi
                        fi

        		if [ $Index_Type = 14 ] ; then
               			echo "    The actual grid size is : ${grid_size}x${grid_size}x${grid_size}"
                		read -p "    The new grid size is : " grid_size
        		fi

        		if [ $Index_Type = 15 ] ; then
                		whiling=''
                		while [[ ! $whiling =~ ^(n|y)$ ]]
                		do
                       	 		read -p "    Automatic save of the generated CUBE-file (y/n) : " whiling
                		done
                		if [ $whiling = n ] ; then
                        		save_cube='No'
                		elif [ $whiling = y ] ; then
                        		save_cube='Yes'
                		fi
        		fi

			if [ $Index_Type = 16 ] ; then
                		echo "    The cutoff for the transition decomposition is : ${cutoff}"
                		read -p "    Please enter the new cutoff : " cutoff
        		fi

			if [ $Index_Type = 17 ] ; then
                                echo "    The name for the Tozer file is by default : Tozer-state"
                                read -p "    Please enter the new suffixe : " name
                        fi
       			
		        if [ $Index_Type = 5 ]; then
				echo "    =============================================================="
                        	echo "         Generation options:"
                        	echo ""
                        	echo "               - From transitions: a"
                        	echo "               - From a list of orbitals: b"
                        	echo "               - From an interval of orbitals: c"   
 				echo ""
				echo "               -Exit: 21"
                        	echo "    --------------------------------------------------------------"
                        	read -p "    Enter a type: " Index_Type
                        	echo "    =============================================================="
                        	if not_in_list "$index_label"  " "  $Index_Type ; then
                                	echo ""
                                	echo "#########################################################################"
                                	echo ""
                                	echo "                 ERROR: Please choose a valid option "
                                	echo ""
                                	echo "#########################################################################"
                                	echo ""
                                
                        	fi
				
				if [ $Index_Type = a ]; then
					
                       	        	read -p "    Excited State Number: " State_Index
                                	echo "    =============================================================="
                                	echo ""	
                                	sed -i 's/nroots/nroots/Ig' ${input_file}
					nb_state=(`python -c"import excited_index; excited_index.find_nb_state_calc('${input_file}')"`)
                                	if [ ${nb_state[1]} = True ] ; then
                                        	max_index=$((${nb_state[0]} + ${nb_state[0]}))
                                	else
                                        	max_index=${nb_state[0]}
                                	fi
                                	if [ ${State_Index} -gt $max_index ] ; then
                                        	echo "----------- ERROR -------------"
                                        	echo ""
                                        	echo "    There is no State ${State_Index} in the file ${input_file}"
                                        	echo ""
                                        	echo "--------------------------------"
                                        	echo ""
                                        
                                	fi

                                	if [ ${State_Index} -gt ${nb_state[0]} ] ; then
                                        	State_Index_temp=$((${State_Index} - ${nb_state[0]}))
                                        	i=1
                                        	temp=`grep -A $i "STATE[[:space:]]*${State_Index_temp}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                                        	echo $temp
                                        	while [[ "$temp" != "" ]];
                                        	do
                                                	echo $temp >> ${transition_file}
                                                	((i++))
                                                	temp=`grep -A $i "STATE[[:space:]]*${State_Index_temp}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                                        	done
                                	else
                                        	i=1
                                        	temp=`grep -A $i "STATE[[:space:]]*${State_Index}:  E=" ${input_file} |tail -n 2 |tail -n 1 `
                                        	echo $temp
                                        	while [[ "$temp" != "" ]];
                                        	do
                                                	echo $temp >> ${transition_file}
                                                	((i++))
                                                	temp=`grep -A $i "STATE[[:space:]]*${State_Index}:  E=" ${input_file} |tail -n 2    |tail -n 1`
                                        	done
                                	fi
					MO_virt=($(python -c"import excited_index; excited_index.molecular_decomposition('virtual', True, 'transition.txt')" | tr -d '[],'))
                                	MO_occ=($(python -c"import excited_index; excited_index.molecular_decomposition('occupied', True, 'transition.txt')" | tr -d '[],'))
                              
                                	version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                	if [ ${version} = 5 ]; then
                                        	for occ in ${MO_occ[@]}; do
                                                	/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF 
                                                	4
                                                	${grid_size}
                                                	2
                                                	$occ
                                                	5
                                                	7
                                                	10
                                                	11
EOF
                                        	done

                                        	for virt in ${MO_virt[@]}; do
                                                	/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                               		4
                                                	${grid_size}
                                                	2
                                                	$virt
                                                	5
                                                	7
                                                	10
                                                	11
EOF
                                        	done

                                	elif [ ${version} = 6 ]; then
                                        	for occ in ${MO_occ[@]}; do
                                                	/Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF 
                                                	4
                                                	${grid_size}
                                                	2
                                                	$occ
                                                	5
                                                	7
                                                	11
                                                	12
EOF
                                        	done

                                        	for virt in ${MO_virt[@]}; do
                                                	/Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                	4
                                                	${grid_size}
                                                	2
                                                	$virt
                                                	5
                                                	7
                                                	11
                                                	12
EOF
                                        	done
                        	        fi
				fi

				if [ $Index_Type = b ]; then
					read -p "    List of orbitals: " List_orb
                                        echo "    =============================================================="
					for orb in ${List_orb}; do
						version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
						if [ ${version} = 5 ]; then
							/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                                	4
                                                	${grid_size}
                                                	2
                                                	$orb
                                                	5
                                                	7
                                                	10
                                                	11
EOF
						elif [ ${version} = 6 ]; then
							/Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                	4
                                                	${grid_size}
                                                	2
                                                	$occ
                                                	5
                                                	7
                                                	11
                                                	12
EOF
						fi
					done

				fi

				if [ $Index_Type = c ]; then
                                        read -p " Starting orbital: " begin_orb
					read -p " Ending orbital: " ending_orb
                                        echo "    =============================================================="
					for (( orb=${begin_orb}; orb<(${ending_orb} + 1); orb ++ )); do
                                                version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                                if [ ${version} = 5 ]; then
                                                        /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                                        4
                                                        ${grid_size}
                                                        2
                                                        $orb
                                                        5
                                                        7
                                                        10
                                                        11
EOF
                                                elif [ ${version} = 6 ]; then
                                                        /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                        4
                                                        ${grid_size}
                                                        2
                                                        $occ
                                                        5
                                                        7
                                                        11
                                                        12
EOF
                                                fi
                                        done

                                fi
			fi

			if [ $Index_Type = 1 ]; then
                		
                		read -p "    Excited State Number: " State_Index
                		echo "    =============================================================="
                		echo ""
				sed -i 's/nroots/nroots/Ig' ${input_file}
                		nb_state=(`python -c"import excited_index; excited_index.find_nb_state_calc('${input_file}')"`)
                		if [ ${nb_state[1]} = True ] ; then
                        		max_index=$((${nb_state[0]} + ${nb_state[0]}))
                		else
                        		max_index=${nb_state[0]}
                		fi
                		if [ ${State_Index} -gt $max_index ] ; then
                        		echo "----------- ERROR -------------"
                        		echo ""
                        		echo "    There is no State ${State_Index} in the file ${input_file}"
                        		echo ""
                        		echo "--------------------------------"
                       			echo ""
                        	
                		fi

                		if [ ${State_Index} -gt ${nb_state[0]} ] ; then
                        		State_Index_temp=$((${State_Index} - ${nb_state[0]}))
                        		i=1
                        		temp=`grep -A $i "STATE[[:space:]]*${State_Index_temp}:  E=" ${input_file} |tail -n 2 |tail -n 1`
					echo $temp
                        		while [[ "$temp" != "" ]];
                        		do
                                		echo $temp >> ${transition_file}
                                		((i++))
                                		temp=`grep -A $i "STATE[[:space:]]*${State_Index_temp}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                        		done
                		else
                        		i=1
                        		temp=`grep -A $i "STATE[[:space:]]*${State_Index}:  E=" ${input_file} |tail -n 2 |tail -n 1 `
                        		echo $temp
					while [[ "$temp" != "" ]];
                        		do
                                		echo $temp >> ${transition_file}
                                		((i++))
                                		temp=`grep -A $i "STATE[[:space:]]*${State_Index}:  E=" ${input_file} |tail -n 2    |tail -n 1`
                        		done
                		fi

                		MO_virt=($(python -c"import excited_index; excited_index.molecular_decomposition('virtual', True, 'transition.txt')" | tr -d '[],'))
                		MO_occ=($(python -c"import excited_index; excited_index.molecular_decomposition('occupied', True, 'transition.txt')" | tr -d '[],'))

				version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
				if [ ${version} = 5 ]; then
                			for occ in ${MO_occ[@]}; do
                        			/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF 
                        			4
                        			${grid_size}
                        			2
                        			$occ
                        			5
                        			7
                        			10
                        			11
EOF
                			done

                			for virt in ${MO_virt[@]}; do
                        			/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                        			4
                        			${grid_size}
                        			2
                        			$virt
                        			5
                        			7
                        			10
                        			11
EOF
                			done

				elif [ ${version} = 6 ]; then
                                        for occ in ${MO_occ[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF 
                                                4
                                                ${grid_size}
                                                2
                                                $occ
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

                                        for virt in ${MO_virt[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                2
                                                $virt
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

				fi

                		python excited_index.py tozer ${filename} ${cutoff} Orca
                		rm  square*
                		rm  mult*
				rm transition.txt

                		if [ ${save_cube} = Yes ] ; then
                        		mkdir -p molecular_orbital
                        		mv *.cube molecular_orbital
                		else
                        		rm *mo*.cube
                		fi

                		mv Tozer.txt Tozer_state_Orca_${State_Index}.txt
				if [ -n "$name" ]; then
                                	mv Tozer_state_Orca_${State_Index}.txt Tozer_state_Orca_${State_Index}_${name}.txt
                                fi

                		rm -r __pycache__

			elif [ $Index_Type = 2 ]; then
                		echo "     Note: To be computed, this index needs two states: an initial state and a fianl state."
                		echo "            If the ground state is the initaial state, please answer 0 for the initial state value." 
                		echo "            Otherwise, write the state number."
                		echo "            This computation may take some time."
                		read -p "    Initial State : " Initial_State
                		read -p "    Final State : " Final_State
                		echo "    =============================================================="
                		echo "" 
                		Dct_file=Dct.txt
                		if [ -f ${Dct_file} ] ; then
                        		echo -n "" > ${Dct_file}
                		fi

                		if [ $Initial_State -ge ${Final_State} ] ; then
                        		echo "       ------- ERROR ---------"
                        		echo ""
                        		echo "          It is impossible to compute the Le Bahers index on this transition"
                        		echo ""
                        		echo "        ----------------------"
                        		echo ""
                        	
                		fi
				sed -i 's/nroots/nroots/Ig' ${input_file}
                		nb_state=(`python -c"import excited_index; excited_index.find_nb_state_calc('${input_file}')"`)
                			
				if [ ${nb_state[1]} = True ] ; then
                        		max_index=$((${nb_state[0]} + ${nb_state[0]}))
                		else
                        		max_index=${nb_state[0]}
                		fi
                		if [ ${Initial_State} -gt $max_index ] ; then
                        		echo "        ------- ERROR ---------"
                        		echo ""
                        		echo "          There is no State ${Initial_State} in the file ${input_file}"
                        		echo "" 
                        		echo "        ----------------------"
                        		echo ""
                        	
                			fi
                		if [ ${Final_State} -gt $max_index ] ; then
                        		echo "        ------- ERROR ---------"
                        		echo ""
                        		echo "          There is no State ${Final_State} in the file ${input_file}"
                        		echo ""
                        		echo "        ----------------------"
                        		echo ""
                        
                		fi

                		if [ $Initial_State = 0 ] ; then
					version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                	
					if [ ${version} = 5 ]; then
                        			/home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                		4
                                		${grid_size}
                                		5
                                		7
                                		6
                                		n
                                		${cis_file}
                                		${Final_State}
                                		12      
EOF
				
					elif [ ${version} = 6 ]; then
						/Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                5
                                                7
                                                6
                                                n
                                                ${cis_file}
                                                ${Final_State}
                                                12
EOF	
					fi

                        		if [ $Final_State -lt 10 ] ; then
                                		mv ${cis_file}dp0${Final_State}.cube diff.cube
                        		else
                                		mv ${cis_file}dp${Final_State}.cube diff.cube
                        		fi
                        		echo""
                        		echo "    =============================================================="
                        		python excited_index.py dct ${filename} 0 ${Final_State}
                        		echo "    =============================================================="

				else

					version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                        if [ ${version} = 5 ]; then
                                                /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                5
                                                7
                                                6
                                                n
                                                ${cis_file}
                                                ${Initial_State} ${Final_State}
                                                12      
EOF

                                        elif [ ${version} = 6 ]; then
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                5
                                                7
                                                6
                                                n
                                                ${cis_file}
                                                ${Initial_State} ${Final_State}
                                                12
EOF
       					fi                        

                        		echo ""
                      			echo "    =============================================================="
                        		python excited_index.py dct ${filename} ${Initial_State} ${Final_State}
                        		echo "    =============================================================="
                		fi
                
				rm square*
				if [ ${save_cube} = Yes ] ; then
                        		mkdir -p electron_density
                        		mv *.cube electron_density
                		else
                        		rm *.cube
                		fi
                		mv Dct.txt Dct_State_${Initial_State}_to_${Final_State}_Orca.txt

				if [ -n $name]; then
                                        mv Dct.txt Dct_State_${Initial_State}_to_${Final_State}_Orca.txt Dct.txt Dct_State_${Initial_State}_to_${Final_State}_${name}_Orca.txt
                                fi
                		rm -r __pycache__

			elif [ ${Index_Type} = 3 ] ; then
                		echo "     Note: To be computed, this index needs two states: an initial singlet state and a final triplet state."
                		echo "            Please, write the state number."
                		echo "            The state number is the quantum indexation (i.e. S1, T1, S2, T2, etc) not the orca compation one."
                		echo "            This computation may take some time."
                		read -p "    Initial State (Singlet): " Initial_State
                		read -p "    Final State (Triplet): " Final_State
                		echo "    =============================================================="
                		echo ""
                		if [ ${Initial_State} = 0 ] ; then
                        		echo " ------- ERROR ---------"
                        		echo ""
                        		echo " Please Chose a Valid State !!"
                        		echo ""
                        		echo " ----------------------"
                        		echo ""
                		fi
                		if [ ${Final_State} = 0 ] ; then
                        		echo " ------- ERROR ---------"
                        		echo ""
                        		echo " Please Chose a Valid State !!"
                        		echo ""
                        		echo " ----------------------"
                        		echo ""
                		fi
				sed -i 's/nroots/nroots/Ig' ${input_file}
                		nb_state=(`python -c"import excited_index; excited_index.find_nb_state_calc('${input_file}')"`)
                		if [ ${nb_state[1]} = True ] ; then
                        		max_index=$((${nb_state[0]} + ${nb_state[0]}))
                		else
                        		max_index=${nb_state[0]}
                		fi
                		if [ ${Initial_State} -gt ${nb_state[0]} ] ; then
                        		echo " ------- ERROR ---------"
                        		echo ""
                        		echo "There is no Singlet State ${Initial_State} in the file ${input_file}"
                        		echo ""
                        		echo " ----------------------"
                        		echo ""
                        	
                		fi
                		echo "    The Singlet State Considered is the S${Initial_State}"
                		i=1
                		temp=`grep -A $i "STATE[[:space:]]*${Initial_State}:  E=" ${input_file} |tail -n 2 |tail -n 1 `
                		while [[ "$temp" != "" ]];
                		do
                        		echo $temp >> ${transition_file}
                        		((i++))
                        		temp=`grep -A $i "STATE[[:space:]]*${Initial_State}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                		done
                		MO_virt_s=($(python -c"import excited_index; excited_index.molecular_decomposition('virtual', True, 'transition.txt')" | tr -d '[],'))
                		MO_occ_s=($(python -c"import excited_index; excited_index.molecular_decomposition('occupied', True, 'transition.txt')" | tr -d '[],'))
			
                                version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                if [ ${version} = 5 ]; then
                                        for occ in ${MO_occ[@]}; do
                                                /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF 
                                                4
                                                ${grid_size}
                                                2
                                                $occ
                                                5
                                                7
                                                10
                                                11
EOF
                                        done

                                        for virt in ${MO_virt[@]}; do
                                                /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                2
                                                $virt
                                                5
                                                7
                                                10
                                                11
EOF
                                        done

                                elif [ ${version} = 6 ]; then
                                        for occ in ${MO_occ[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF 
                                                4
                                                ${grid_size}
                                                2
                                                $occ
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

                                        for virt in ${MO_virt[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                2
                                                $virt
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

                                fi		

                		mv transition.txt transition_S${Initial_State}.txt
                		transition_file=transition.txt
                		if [ -f ${transition_file} ] ; then
                        		echo -n "" > ${transition_file}
                		fi
                		if [ ${Final_State} -lt ${nb_state[0]} ] ; then
                        		final_state=${Final_State}
                        		echo ""
                        		echo "    The Triplet State Considered is the T${final_State}"
                        		echo ""
                		elif [ ${Final_State} -gt ${nb_state[0]} ] ; then
                        		if [ ${Final_State} -le ${max_index} ] ; then
                                		final_state=$((${Final_State} - ${nb_state[0]}))
                                		echo "    The Triplet State Considered is th T${final_state}"
                                		echo ""
                        			else
                                		echo " ------- ERROR ---------"
                                		echo ""
                                		echo "There is no Triplet State ${final_state} in the file ${input_file}"
                                		echo ""
                                		echo " ----------------------"
                        		fi
                		fi
                		i=1
                		temp=`grep -A $i "STATE[[:space:]]*${final_state}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                		while [[ "$temp" != "" ]] ;
                		do
                        		echo $temp >> ${transition_file}
                        		((i++))
                        		temp=`grep -A $i "STATE[[:space:]]*${final_state}:  E=" ${input_file} |tail -n 2 |tail -n 1`
                		done
                		MO_virt_t=($(python -c"import excited_index; excited_index.molecular_decomposition('virtual', True, 'transition.txt')" | tr -d '[],'))
                		MO_occ_t=($(python -c"import excited_index; excited_index.molecular_decomposition('occupied', True, 'transition.txt')" | tr -d '[],'))

                                version=$(grep 'Program Version' ${input_file} | sed -n 's/.*Program Version\s*\([0-9]\).*/\1/p')
                                if [ ${version} = 5 ]; then
                                        for occ in ${MO_occ[@]}; do
                                                /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF 
                                                4
                                                ${grid_size}
                                                2
                                                $occ
                                                5
                                                7
                                                10
                                                11
EOF
                                        done

                                        for virt in ${MO_virt[@]}; do
                                                /home/ssteinma/softs/orca_5_0_4_linux_x86-64_openmpi411/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                2
                                                $virt
                                                5
                                                7
                                                10
                                                11
EOF
                                        done

                                elif [ ${version} = 6 ]; then
                                        for occ in ${MO_occ[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF 
                                                4
                                                ${grid_size}
                                                2
                                                $occ
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

                                        for virt in ${MO_virt[@]}; do
                                                /Xnfs/chimie/debian11/orca/orca_6_0_1/orca_plot ${gbw_file} -i << EOF
                                                4
                                                ${grid_size}
                                                2
                                                $virt
                                                5
                                                7
                                                11
                                                12
EOF
                                        done

                                fi

                		mv transition.txt transition_T${final_state}.txt
                		echo ""
                		echo "    =============================================================="
                		python excited_index.py omega_soc ${filename} ${cutoff} ${Initial_State} ${final_state}
                		echo "    =============================================================="
                		rm square*
                		rm mult*
                		rm transition*
                
				if [ ${save_cube} = Yes ] ; then
                        		mkdir -p molecular_orbital
                        		mv *.cube molecular_orbital
                		else
                        		rm *.cube
                		fi
                		mv Omega.txt Omega_S${Initial_State}_T${Final_State}_Orca.txt

				if [ -n $name]; then
                                        mv Omega.txt Omega_S${Initial_State}_T${Final_State}_Orca.txt Omega.txt Omega_S${Initial_State}_T${Final_State}_${name}_Orca.txt
                                fi
                		rm -r __pycache__

			elif [ $Index_Type = 4 ]; then
                		echo "     Note : This index is still in implementation"
                		echo "            Please give the CUBE-file of the electron density of the final and the initial state ."
                		read -p "    CUBE-file of the Electron Density of the Initial State : " Initial_State_file
                		read -p "    CUBE-file of the Electron Density of the Final State : " Final_State_file
                		echo "    =============================================================="
                		python excited_index.py double_overlap_cube ${Initial_State_file} ${Final_State_file}
                		echo "    =============================================================="
                		rm mult.cube
                		rm -rf __pycache__
        		fi
		done
	fi	
done
echo  "                         ,"
echo  "  ,-.       _,---._ __  / \\"
echo  " /  )    .-'       \`./ /   \\"
echo  "(  (   ,'            \`/    /|"
echo  " \\  \`-\"             \\'\\   / |"
echo  "  \`.              ,  \\ \\ /  |"
echo  "   /\`.          ,'-\`----Y   |"
echo  "  (            ;        |   '"
echo  "  |  ,-.    ,-'         |  /"
echo  "  |  | (   |            | /"
echo  "  )  |  \\  \`.___________|/"
echo  "  \`--'   \`--'"
echo ""
echo "                      Normal Terminaison"
echo ""
echo "    =============================================================="

