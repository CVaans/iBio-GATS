#!/bin/sh
echo "A workflow for insect odorant receptor structural modelling..."

#Change directory to the ibiogats root directory
cd ./ibiogats

#Read the input target sequence in fasta file format
#The output showing the score, resolution, sequence identity for each best template selected will be displayed here 
#The result summary for each best template will be available for download for user reference
python3.7 ./ibiogats/main.py sequence.fasta

selectedTemplate="`cat ./ibiogats/selected-template.txt`"
echo "Template selected for this sequence... $selectedTemplate"

# Change directory for Modeller files
cd ./modeller-files/$selectedTemplate

receptorArray=()

#Check for modeller input files
for f in *; do
    if [ -d "$f" ]; then
        # Will not run if no directories are available
        receptorArray[ ${#receptorArray[@]} ]=$f
    fi
done

#Get the user to select any one of the 3 structures (7LIG,7LID,7LIC) if MhOR5 is selected for model building
#else strcuture (6C70) is selected if user selects Orco for model building
if [ ${#receptorArray[@]} -gt 0 ]; then
    count=1
    for model in ${receptorArray[@]}; do
        echo "$count. $model"
        count=$(( $count + 1 ))
    done

    echo "Please choose any one of the three structure for the template MhOR5:"
    read userOption
    userOption=$(( $userOption - 1 ))
    bestTemplate=${receptorArray[$userOption]}
    cd ./modeller-files/$selectedTemplate/$bestTemplate
    python3.7 ./modeller-files/$selectedTemplate/$bestTemplate/model.py
else
    python3.7 ./modeller-files/$selectedTemplate/model.py 
fi


