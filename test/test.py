import unittest
import capsid
import logging
from mock import patch

class Object(object):
    pass

class MockConfig(object):

    configuration = Object()
    configuration.username = ""
    configuration.password = ""
    configuration.host = "localhost"
    configuration.port = "27017"
    configuration.database = "test_capsid"
 
    def get(self, key1, key2):
        return getattr(self.configuration, key2)

class CapsidTest(unittest.TestCase):

    args = Object()
    db = None

    def setUp(self):
        logging.basicConfig(level=logging.INFO)
        self.args.logging = logging

        configuration = MockConfig()
        setattr(capsid.database, 'get_configuration', lambda: configuration)

        self.db = capsid.database.connect(self.args)
        capsid.configure.logger = logging
        capsid.configure.db = self.db

        self.db = None


    def tearDown(self):
        if self.db:
            self.db.project.remove()
            self.db.role.remove()
            self.db.connection.disconnect()


    def test_database_connection(self):
        '''
        Check that a database connection is truly possible
        '''
        self.db = capsid.database.connect(self.args)
        assert self.db


    def test_ensure_indexes(self):
        '''
        Check that we can create all the right indexes. This should be
        done fairly early on. 
        '''

        logger = Object()
        logger.debug = lambda string: string

        self.db = capsid.database.connect(self.args)
        capsid.configure.logger = logger
        capsid.configure.db = self.db
        capsid.configure.ensure_indexes()

        # This doesn't return a meaningful value, but we ought to be able to find
        # some evidence of an index. To do.
        pass


    def test_create_project(self):
        '''
        Check that we can create a project
        '''

        self.db = capsid.database.connect(self.args)

        setattr(self.args, 'pdesc', "Project description")
        setattr(self.args, 'project', "morag")
        setattr(self.args, 'pname', "Project name")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)
        
        # We now use mongo to check we have that project
        project = self.db.project.find_one({"label": "morag"})
        assert project
        assert project['description'] == "Project description"
        assert project['label'] == "morag"
        assert project['name'] == "Project name"
        assert project['wikiLink'] == "Project link"


    def test_create_duplicate_project(self):
        '''
        Check that we can create a project
        '''

        self.db = capsid.database.connect(self.args)

        setattr(self.args, 'pdesc', "First description")
        setattr(self.args, 'project', "morag")
        setattr(self.args, 'pname', "First name")
        setattr(self.args, 'link', "First link")
        capsid.project.main(self.args)
        
        # Now make a second one, with the name project key but different values
        setattr(self.args, 'pdesc', "Second description")
        setattr(self.args, 'project', "morag")
        setattr(self.args, 'pname', "Second name")
        setattr(self.args, 'link', "Second link")

        try:
            capsid.project.main(self.args)
        except SystemExit as inst:
            assert inst.code == 1, "Should have an exit status of 1"


if __name__ == '__main__':
    unittest.main()