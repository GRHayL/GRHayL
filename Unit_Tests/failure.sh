#!/bin/sh

if $1; then
    echo "Failed to fail!"
    exit 1
else
    echo "Failed successfully!"
fi
