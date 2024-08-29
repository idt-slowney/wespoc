#!/bin/bash

# Prompt
echo "Thanks for using the IDT Nextflow Workflow Template!"
echo "Your workflow has been initialized with default values."
read -p "Please provide a name for your workflow: " replacement

if [ -z "$replacement" ]; then
    echo "Please enter a valid workflow name."
    exit 1
fi

replacement_upper=$(echo "$replacement" | tr '[:lower:]' '[:upper:]')

find . -type f -exec sed -i \
    -e "s/wespoc/$replacement/gI" \
    -e "s/wespoc/$replacement/g" {} +

find . -depth -name "*wespoc*" | while read -r filename; do 
    new_name=$(echo "$filename" | sed \
        -e "s/wespoc/$replacement/gI" \
        -e "s/wespoc/$replacement/g")
    mv "$filename" "$new_name"
done

echo "Your workflow directory has been updated with your workflow name."