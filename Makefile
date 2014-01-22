ADD=README.md Makefile
COM=first commit

gadd:
	git add $(ADD)
gcom:
	git commit -m "$(COM)"
gcomam:
	git commit --amend
gstat:
	git status
gpush:
	git push -u origin master

