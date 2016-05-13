#!/bin/bash

IFS=','
while read -r name repo tag path commands || [[ -n $name ]]; do
	name=`echo $name | xargs`
	repo=`echo $repo | xargs`
	tag=`echo $tag | xargs`
	path=`echo $path | xargs`
	#echo $name
	#echo $repo
	#echo $tag
	#echo $path
	#echo $commands
	if [ -f $name.version ]; then
		if [ "`cat $name.version`x" = "${tag}x" ]; then
			continue
		fi
	fi
	if [ -d $name ]; then
		rm -rf $name
	fi
	mkdir $name
	rm -rf .tmp
	git clone $repo .tmp
	cd .tmp
	git reset --hard $tag
	rm -rf .git
	cd ..
	cp -R .tmp$path/* $name
	rm -rf .tmp
	cd $name
	for cmd in $commands; do
		cmd=`echo $cmd | xargs`
		#echo $cmd
		eval `echo $cmd`
	done
	cd ..
	echo $tag > $name.version
done < externals
