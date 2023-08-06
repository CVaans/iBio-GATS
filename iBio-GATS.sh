#!/bin/sh
echo "**************************************************************"
echo "iBio-GATS, a semi-automated tool for insect OR structural modelling" 
echo "**************************************************************"

echo "Checking and processing the given input sequence"

#Change directory to the ibiogats root directory
cd ./ibiogats

#Read the input target sequence in fasta file format
#The output showing the template name, SSD score, resolution for each best template selected will be displayed here 
#The result summary for each best template will be available for download for user reference
python3.7 ./ibiogats/main.py sequence.fasta

selectedTemplate="`cat ./ibiogats/selected-template.txt`"
echo "Template selected for this sequence... $selectedTemplate"

#Copy selected template to modeller folder
cp ./ibiogats/$selectedTemplate.ali ./modeller-files/

template_text=$selectedTemplate

# Set the delimiter
IFS="-"

# Split the text and store it in an array
read -ra file_names <<< "$template_text"
$./ibiogats/deactivate

# folder path for Modeller files
IFS=""
cd ./modeller-files/

# user selected template and its related files were sent to Modeller 10.4 
    
ali_file_name=$selectedTemplate.ali
pdb_file_name=${file_names[0]}
target_ali_file_name=${file_names[1]}

#printing the name of selected template and its related files for user reference

echo "The following files are used for model building"
echo "Target file name: $target_ali_file_name"
echo "PDB file name: $pdb_file_name"
echo "Template Target alignment file name: $ali_file_name"
echo "Please check the alignment file. Do you wish to manually edit the alignment file"

# prompting the user for using the option of manually editing the template-target alignment file.

while true; do
    read -p "Please answer Y or N   " yn
    case $yn in
        [Yy]* ) echo "Ok. Please change the alignment and save the file";
                 read -p "Please press any key and enter to proceed after you finish editing the alignment file.      ";
                 echo "Proceeding to model building ";
                 sleep 5s;
                 break;;
        [Nn]* ) echo "Continuing with the same alignment ";
                echo "Proceeding to model building ";
                break;;
        * ) echo "Invalid response";;
    esac
    
done
sleep 10s
#running model.py in Modeller 10.4 for target file with the selected template and the template-target alignment file
python3.7 /mnt/c/vaanathi/Ubuntu-Bio-Gats/modeller-files/model.py $ali_file_name $pdb_file_name $target_ali_file_name
pwd
echo "Model building completed for $target_ali_file_name using the template $pdb_file_name"
echo "Please save the above models in seperate folders for your reference"