ADD=README.md Makefile

gadd:
	git add $(ADD)
gcom:
	git commit -m "first commit"
gcomam:
	git commit --amend
gstat:
	git status
gpush:
	git push -u origin master

