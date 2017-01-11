# MB-DB: MetaBolomic DataBase

## Installation

These instructions are written for Ubuntu 16.04. If you require any platform-specific assistance in setting MB-DB up, please feel free to [open an issue](https://github.com/KeironO/MB-DB/issues).

### MongoDB

#### Installing MongoDB

MongoDB is already included in the standard Ubuntu package repository, but the recommended installation source is from the official MongoDB repository.

To do this, firstly you are required to import the key for the official MongoDB repository to your sources list.

```
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv EA312927
```

Now you are required to add the repository details through the addition of the packages list.

```
echo "deb http://repo.mongodb.org/apt/ubuntu xenial/mongodb-org/3.2 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-3.2.list && sudo apt-get update
```
Once completed, you can now install the MongoDB package.

```
sudo apt-get install -y mongodb-org
```

#### Configurating MongoDB

We're not quite finished(!) We now need to tell ```systemd``` how to manage the resource.

First of all, create a new ```service``` file in the ```/etc/systemd/system/``` directory.

```
sudo vim /etc/systemd/system/mongodb.service
```

Paste in the following text.

```
[Unit]
Description=High-performance, schema-free document-oriented database
After=network.target

[Service]
User=mongodb
ExecStart=/usr/bin/mongod --quiet --config /etc/mongod.conf

[Install]
WantedBy=multi-user.target
```
Then save and close the file (for the vim virgins to do this press ```Esc``` and then type ```:wq```).

Now set up the service to run by default.

```
sudo systemctl start mongodb && sudo systemctl enable mongodb
```

#### Setting up the collection

```
mongoimport --db mbdb --collection metabolites --type json --file mb-db.json --jsonArray 
```

### python-pip

## License

Code released under [the MIT license](https://github.com/MB-DB/omicsdb/blob/master/LICENSE).
