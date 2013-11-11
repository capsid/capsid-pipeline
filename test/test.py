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
        logging.basicConfig(level=logging.FATAL)
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
        self.assertIsNotNone(self.db)


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
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'pname', "Simulated")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)
        
        # We now use mongo to check we have that project
        project = self.db.project.find_one({"label": "simu"})

        self.assertIsNotNone(project)
        self.assertEqual(project['description'], "Project description")
        self.assertEqual(project['label'], "simu")
        self.assertEqual(project['name'], "Simulated")
        self.assertEqual(project['wikiLink'], "Project link")


class CapsidProjectTest(unittest.TestCase):

    args = Object()
    db = None

    def setUp(self):
        logging.basicConfig(level=logging.FATAL)
        self.args.logging = logging

        configuration = MockConfig()
        setattr(capsid.database, 'get_configuration', lambda: configuration)

        self.db = capsid.database.connect(self.args)
        capsid.configure.logger = logging
        capsid.configure.db = self.db

        self.db = capsid.database.connect(self.args)

        setattr(self.args, 'pdesc', "Project description")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'pname', "Simulated")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)

        setattr(self.args, 'pdesc', "Other project description")
        setattr(self.args, 'project', "other")
        setattr(self.args, 'pname', "Other")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)


    def tearDown(self):
        self.db.project.remove()
        self.db.role.remove()
        self.db.sample.remove()
        self.db.connection.disconnect()


    def test_create_duplicate_project(self):
        '''
        Check that we can't create a second project with the same label
        '''

        # Now make a second one, with the name project key but different values
        setattr(self.args, 'pdesc', "Second description")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'pname', "Second name")
        setattr(self.args, 'link', "Second link")

        exitCode = 0

        try:
            capsid.project.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 1, "Should have an exit status of 1")


    def test_create_sample(self):
        '''
        Check that we can create a sample.
        '''

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Simulated Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")

        exitCode = 0

        try:
            capsid.sample.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 0, "Should have an exit status of 0")

        sample = self.db.sample.find_one({"name": "SIMU001"})

        self.assertIsNotNone(sample)
        self.assertEqual(sample['name'], "SIMU001")
        self.assertEqual(sample['cancer'], "simulated")
        self.assertEqual(sample['description'], "Simulated Sample")
        self.assertEqual(sample['role'], "CASE")
        self.assertEqual(sample['source'], "P")
        self.assertEqual(sample['projectLabel'], "simu")

        self.assertIsNotNone(sample['projectId'])

        project = self.db.project.find_one({"_id" : sample['projectId']})
        self.assertIsNotNone(project)
        self.assertEqual(project['label'], "simu")


    def test_create_duplicate_sample(self):
        '''
        Check that we can't create a sample with the same identifier, within
        the same project.
        '''

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Simulated Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")

        capsid.sample.main(self.args)

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Second Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")

        exitCode = 0

        try:
            capsid.sample.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 1, "Should have an exit status of 1")


    def test_create_sample_second_project(self):
        '''
        Check that we can create a sample with the same identifier, within
        a different project.
        '''

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Simulated Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")

        capsid.sample.main(self.args)

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "other")
        setattr(self.args, 'sdesc', "Second Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")

        exitCode = 0

        try:
            capsid.sample.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 0, "Should have an exit status of 0")



class CapsidSampleTest(unittest.TestCase):

    args = Object()
    db = None

    def setUp(self):
        logging.basicConfig(level=logging.FATAL)
        self.args.logging = logging

        configuration = MockConfig()
        setattr(capsid.database, 'get_configuration', lambda: configuration)

        self.db = capsid.database.connect(self.args)
        capsid.configure.logger = logging
        capsid.configure.db = self.db

        self.db = capsid.database.connect(self.args)

        setattr(self.args, 'pdesc', "Project description")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'pname', "Simulated")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)

        setattr(self.args, 'pdesc', "Other project description")
        setattr(self.args, 'project', "other")
        setattr(self.args, 'pname', "Other")
        setattr(self.args, 'link', "Project link")
        capsid.project.main(self.args)

        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Simulated Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")
        capsid.sample.main(self.args)

        setattr(self.args, 'sample', "SIMU002")
        setattr(self.args, 'project', "simu")
        setattr(self.args, 'sdesc', "Second Sample")
        setattr(self.args, 'cancer', "simulated")
        setattr(self.args, 'source', "P")
        setattr(self.args, 'role', "CASE")
        capsid.sample.main(self.args)


    def tearDown(self):
        self.db.project.remove()
        self.db.role.remove()
        self.db.sample.remove()
        self.db.alignment.remove()
        self.db.connection.disconnect()


    def test_create_alignment(self):
        '''
        Check that we can create an alignment.
        '''

        setattr(self.args, 'aligner', "Aligner")
        setattr(self.args, 'platform', "Platform")
        setattr(self.args, 'align', "simulated")
        setattr(self.args, 'type', "Type")
        setattr(self.args, 'projectLabel', "simu")
        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'infile', "infile")
        setattr(self.args, 'outfile', "outfile")

        exitCode = 0

        try:
            capsid.alignment.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 0, "Should have an exit status of 0")

        alignment = self.db.alignment.find_one({"name": "simulated"})
        self.assertIsNotNone(alignment)
        self.assertEqual(alignment["aligner"], "Aligner")
        self.assertEqual(alignment["platform"], "Platform")
        self.assertEqual(alignment["name"], "simulated")
        self.assertEqual(alignment["type"], "Type")
        self.assertEqual(alignment["projectLabel"], "simu")
        self.assertEqual(alignment["sample"], "SIMU001")
        self.assertEqual(alignment["infile"], "infile")
        self.assertEqual(alignment["outfile"], "outfile")
        self.assertIsNotNone(alignment["sampleId"])
        self.assertIsNotNone(alignment["projectId"])


    def test_create_duplicate_alignment(self):
        '''
        Check that we can't create an alignment with the same name in the same sample.
        '''

        setattr(self.args, 'aligner', "Aligner")
        setattr(self.args, 'platform', "Platform")
        setattr(self.args, 'align', "simulated")
        setattr(self.args, 'type', "Type")
        setattr(self.args, 'projectLabel', "simu")
        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'infile', "infile")
        setattr(self.args, 'outfile', "outfile")
        capsid.alignment.main(self.args)


        setattr(self.args, 'aligner', "Second Aligner")
        setattr(self.args, 'platform', "Second Platform")
        setattr(self.args, 'align', "simulated")
        setattr(self.args, 'type', "Second Type")
        setattr(self.args, 'projectLabel', "simu")
        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'infile', "Second infile")
        setattr(self.args, 'outfile', "Second outfile")

        exitCode = 0

        try:
            capsid.alignment.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 1, "Should have an exit status of 1")


    def test_create_alignment_second_sample(self):
        '''
        Check that we can create an alignment with the same name, so long as 
        do that in a different sample.
        '''

        setattr(self.args, 'aligner', "Aligner")
        setattr(self.args, 'platform', "Platform")
        setattr(self.args, 'align', "simulated")
        setattr(self.args, 'type', "Type")
        setattr(self.args, 'projectLabel', "simu")
        setattr(self.args, 'sample', "SIMU001")
        setattr(self.args, 'infile', "infile")
        setattr(self.args, 'outfile', "outfile")
        capsid.alignment.main(self.args)


        setattr(self.args, 'aligner', "Second Aligner")
        setattr(self.args, 'platform', "Second Platform")
        setattr(self.args, 'align', "simulated")
        setattr(self.args, 'type', "Second Type")
        setattr(self.args, 'projectLabel', "simu")
        setattr(self.args, 'sample', "SIMU002")
        setattr(self.args, 'infile', "Second infile")
        setattr(self.args, 'outfile', "Second outfile")

        exitCode = 0

        try:
            capsid.alignment.main(self.args)
        except SystemExit as inst:
            exitCode = inst.code

        self.assertEqual(exitCode, 0, "Should have an exit status of 0")


if __name__ == '__main__':
    unittest.main()