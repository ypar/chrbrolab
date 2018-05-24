
# README for lab-wide communications and other setup trivia
# YoSon Park



## Index

<!--ts-->

* [slack](#slack)
* [google drive](#gdrive)
* [google cloud platform](#gcp)
* [dropbox](#dropbox)
* [vpn](#vpn)

<!--te-->


I have setup some docs and non-pmacs workspace for convenience of chrbrolab members. I no longer maintain any of these but I list a few that will continue to be useful for incoming members of the group.



# slack

I have set up a slack workspace for chrbrolab members. If you are new to the group, cdb and mc are now admins and can give you access.




# gdrive

There are several documentations and meeting schedules, etc. maintained and shared in the lab gdrive. Contact cdb to give you access. 



# gcp

I have made a lab account for a shared google cloud platform access. This is useful for collaborative projects. Contact me and cdb if you want to ask for access.




# dropbox

This is seldom used but will be useful if you have lab protocols, etc. that need to be shared. Contact gh to give you the user account information.





# vpn

Here is an instruction for vpn setup for pmacs off campus:

```
SETTING VPN
===========
1. Download and install Forticlient from http://forticlient.com/#download
2. After successful installation, open application as:
  Go —> Applications —> Forticlient
3. Within FortiClient application, setup credentials as:
  (a) Click “Remote Access”
  (b) ONE time only: click settings icon next to [VPN Name] field —> select "Add a new connection" and fill information as:
      VPN Type —> SSL VPN
      Connection Name —> PMACS
      Description —>
      Remote Gateway —> juneau.med.upenn.edu
                      port: 443
      Client Certificate —> None
      Authentication —>
      Username —> <your-pmacs-username>
      Apply
4. Provide just created credentials in the “Remote Access” tab
  VPN Name —> PMACS
  Username —> <your-pmacs-username>
  Password —> <your-pmacs-password>
  Connect
5. Open terminal and ssh to PMACS:
  ssh username@consign.pmacs.upenn.edu
```

