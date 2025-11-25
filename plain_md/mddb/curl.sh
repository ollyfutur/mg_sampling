#!/bin/bash

# Check if no arguments are provided or if the number of arguments is odd
if [ $# -eq 0 ] || [ $(( $# % 2 )) -ne 0 ]; then
  echo "Error: No arguments provided or arguments must be provided in pairs (file path and filename)."
  echo "Usage: $0 <file path1> <filename1> <file path2> <filename2> ..."
  exit 1
fi

echo "Uploading files. Depending on the file(s) size, this may take a while."

# Loop through the pairs of arguments
while [ $# -gt 0 ]; do
  file=$1
  filename=$2

  resource="/20250530/94ehml/${filename}"
  date=`date -R`
  request="PUT\n\napplication/octet-stream\n${date}\n${resource}"
  signature=`echo -en ${request} | openssl sha1 -hmac 5Au7bJ8iHS7PygcHBq769DvWfFa7haZ+MpLpolX7 -binary | base64`

  curl --progress-bar --upload-file "${file}" \
      -H "Host: jsc.mddbr.eu:7000" \
      -H "Date: ${date}" \
      -H "Content-Type: application/octet-stream" \
      -H "Authorization: AWS GZJ11DUZAI8SGTY1PH8O:${signature}" \
      https://jsc.mddbr.eu:7000${resource} | type

  # Shift to the next pair of arguments
  shift 2
done

echo "Done!"
