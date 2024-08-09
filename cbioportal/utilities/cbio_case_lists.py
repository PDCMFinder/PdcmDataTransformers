import sys
from os import makedirs
from os.path import isfile, join, isdir, exists
class cbioportal_case_lists:
    def __init__(self, filename, study_dir, study, overwrite, verbose):
        self.config_file_name = filename
        self.case_list_dir = join(study_dir, "case_lists")
        if not exists(self.case_list_dir):
            makedirs(self.case_list_dir)
        self.study_dir = study_dir
        self.study_id = study
        self.overwrite = overwrite
        self.verbose = verbose
        self.CASE_LIST_CONFIG_HEADER_COLUMNS = ["CASE_LIST_FILENAME", "STAGING_FILENAME", "META_STABLE_ID",
                                                "META_CASE_LIST_CATEGORY", "META_CANCER_STUDY_ID",
                                                "META_CASE_LIST_NAME",
                                                "META_CASE_LIST_DESCRIPTION"]
        self.CASE_LIST_UNION_DELIMITER = "|"
        self.CASE_LIST_INTERSECTION_DELIMITER = "&"
        self.MUTATION_STAGING_GENERAL_PREFIX = "data_mutations"
        self.SEQUENCED_SAMPLES_FILENAME = "sequenced_samples.txt"
        self.MUTATION_CASE_LIST_META_HEADER = "sequenced_samples"
        self.MUTATION_CASE_ID_COLUMN_HEADER = "Tumor_Sample_Barcode"
        self.SAMPLE_ID_COLUMN_HEADER = "SAMPLE_ID"
        self.NON_CASE_IDS = frozenset(
            ["MIRNA", "LOCUS", "ID", "GENE SYMBOL", "ENTREZ_GENE_ID", "HUGO_SYMBOL", "LOCUS ID", "CYTOBAND",
             "COMPOSITE.ELEMENT.REF", "HYBRIDIZATION REF"])
        self.CANCER_STUDY_TAG = "<CANCER_STUDY>"
        self.NUM_CASES_TAG = "<NUM_CASES>"

    def generate_case_lists(self):
        header = []
        with open(self.config_file_name, 'r') as case_list_config_file:
            # get header and validate
            header = case_list_config_file.readline().rstrip('\n').rstrip('\r').split('\t')
            # check full header matches what we expect
            for column in self.CASE_LIST_CONFIG_HEADER_COLUMNS:
                if column not in header:
                    print >> sys.stderr, "ERROR: column '%s' is not in '%s'" % (column, self.config_file_name)
                    sys.exit(2)

            for line in case_list_config_file:
                line = line.rstrip('\n').rstrip('\r')
                config_fields = line.split('\t')
                case_list_filename = config_fields[header.index("CASE_LIST_FILENAME")]
                staging_filename_list = config_fields[header.index("STAGING_FILENAME")]
                case_list_file_full_path = join(self.case_list_dir, case_list_filename)
                if isfile(case_list_file_full_path) and not self.overwrite:
                    if self.verbose:
                        print("LOG: generate_case_lists(), '%s' exists and overwrite is false, skipping caselist..." % (
                            case_list_filename))
                    continue

                # might be single staging file
                staging_filenames = []
                # union (like all cases)
                union_case_list = self.CASE_LIST_UNION_DELIMITER in staging_filename_list
                # intersection (like complete or cna-seq)
                intersection_case_list = self.CASE_LIST_INTERSECTION_DELIMITER in staging_filename_list
                delimiter = self.CASE_LIST_UNION_DELIMITER if union_case_list else self.CASE_LIST_INTERSECTION_DELIMITER
                staging_filenames = staging_filename_list.split(delimiter)
                if self.verbose:
                    print("LOG: generate_case_lists(), staging filenames: %s" % (",".join(staging_filenames)))

                # if this is intersection all staging files must exist
                if intersection_case_list and \
                        not all([isfile(join(self.study_dir, intersection_filename)) for intersection_filename in
                                 staging_filenames]):
                    continue

                # this is the set we will pass to write_case_list_file
                case_set = set([])
                # this indicates the number of staging files processed -
                # used to verify that an intersection should be written
                num_staging_files_processed = 0
                for staging_filename in staging_filenames:
                    if self.verbose:
                        print("LOG: generate_case_lists(), processing staging file '%s'" % (staging_filename))
                    # compute the case set
                    case_list = []
                    case_list = self.get_case_list_from_staging_file(staging_filename)

                    if len(case_list) == 0:
                        if self.verbose:
                            print("LOG: generate_case_lists(), no cases in '%s', skipping..." % (staging_filename))
                        continue

                    if intersection_case_list:
                        if len(case_set) == 0:
                            # it is empty so initialize it
                            case_set = set(case_list)
                        else:
                            case_set = case_set.intersection(case_list)
                    else:
                        # union of files or single file
                        case_set = case_set.union(case_list)

                    num_staging_files_processed += 1

                # write case list file (don't make empty case lists)
                if len(case_set) > 0:
                    if self.verbose:
                        print("LOG: generate_case_lists(), calling write_case_list_file()...")

                    # do not write out complete cases file unless we've processed all the files required
                    if intersection_case_list and num_staging_files_processed != len(staging_filenames):
                        if self.verbose:
                            print(
                                "LOG: generate_case_lists(), number of staging files processed (%d) != number of staging files required (%d) for '%s', skipping call to write_case_list_file()..." % (
                                    num_staging_files_processed, len(staging_filenames), case_list_filename))
                    else:
                        self.write_case_list_file(header, config_fields, case_list_file_full_path, case_set)
                elif self.verbose:
                    print(
                        "LOG: generate_case_lists(), case_set.size() == 0, skipping call to write_case_list_file()...")

    def get_case_list_from_staging_file(self, staging_filename):
        if self.verbose:
            print("LOG: get_case_list_from_staging_file(), '%s'" % (staging_filename))

        case_set = set([])

        # if we are processing mutations data and a SEQUENCED_SAMPLES_FILENAME exists, use it
        if self.MUTATION_STAGING_GENERAL_PREFIX in staging_filename:
            sequenced_samples_full_path = join(self.study_dir, self.SEQUENCED_SAMPLES_FILENAME)
            if isfile(sequenced_samples_full_path):
                if self.verbose:
                    print(
                        "LOG: get_case_list_from_staging_file(), '%s' exists, calling get_case_list_from_sequenced_samples_file()" % (
                            self.SEQUENCED_SAMPLES_FILENAME))
                return self.get_case_list_from_sequenced_samples_file(sequenced_samples_full_path)

        staging_file_full_path = join(self.study_dir, staging_filename)
        if not isfile(staging_file_full_path):
            return []

        # staging file
        with open(staging_file_full_path, 'r') as staging_file:
            id_column_index = 0
            process_header = True
            for line in staging_file:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    if line.startswith('#' + self.MUTATION_CASE_LIST_META_HEADER + ':'):
                        # split will split on any whitespace, tabs or any number of consecutive spaces
                        return line[len(self.MUTATION_CASE_LIST_META_HEADER) + 2:].strip().split()
                    continue  # this is a comment line, skip it
                values = line.split('\t')

                # is this the header line?
                if process_header:
                    # look for MAF file case id column header
                    # if this is not a MAF file and header contains the case ids, return here
                    # we are assuming the header contains the case ids because SAMPLE_ID_COLUMN_HEADER is missing
                    if self.MUTATION_CASE_ID_COLUMN_HEADER not in values and self.SAMPLE_ID_COLUMN_HEADER not in [
                        x.upper() for x
                        in
                        values]:
                        if self.verbose:
                            print(
                                "LOG: get_case_list_from_staging_file(), this is not a MAF header but has no '%s' column, we assume it contains sample ids..." % (
                                    self.SAMPLE_ID_COLUMN_HEADER))
                        for potential_case_id in values:
                            # check to filter out column headers other than sample ids
                            if potential_case_id.upper() in self.NON_CASE_IDS:
                                continue
                            case_set.add(potential_case_id)
                        break  # got case ids from header, don't read the rest of the file
                    else:
                        # we know at this point one of these columns exists, so no fear of ValueError from index method
                        id_column_index = values.index(
                            self.MUTATION_CASE_ID_COLUMN_HEADER) if self.MUTATION_CASE_ID_COLUMN_HEADER in values else [
                            x.upper()
                            for
                            x in
                            values].index(
                            self.SAMPLE_ID_COLUMN_HEADER)
                        if self.verbose:
                            print(
                                "LOG: get_case_list_from_staging_file(), this is a MAF or clinical file, samples ids in column with index: %d" % (
                                    id_column_index))
                    process_header = False
                    continue  # done with header, move on to next line
                case_set.add(values[id_column_index])

        return list(case_set)

    def get_case_list_from_sequenced_samples_file(self, sequenced_samples_full_path):
        if self.verbose:
            print("LOG: get_case_list_from_sequenced_samples_file, '%s'", sequenced_samples_full_path)

        case_set = set([])
        with open(sequenced_samples_full_path, 'r') as sequenced_samples_file:
            for line in sequenced_samples_file:
                case_set.add(line.rstrip('\n'))

        if self.verbose:
            print("LOG: get_case_list_from_sequenced_samples_file, case set size: %d" % (len(case_set)))

        return list(case_set)

    def write_case_list_file(self, case_list_config_header, case_list_config_fields, case_list_full_path, case_set):
        if self.verbose:
            print("LOG: write_case_list_file(), '%s'" % (case_list_full_path))
        with open(case_list_full_path, 'w') as case_list_file:
            case_list_file.write("cancer_study_identifier: " + self.study_id + "\n")
            stable_id = case_list_config_fields[case_list_config_header.index("META_STABLE_ID")].replace(
                self.CANCER_STUDY_TAG,
                self.study_id)
            case_list_file.write("stable_id: " + stable_id + "\n")
            case_list_file.write(
                "case_list_name: " + case_list_config_fields[
                    case_list_config_header.index("META_CASE_LIST_NAME")] + "\n")
            case_list_description = case_list_config_fields[
                case_list_config_header.index("META_CASE_LIST_DESCRIPTION")].replace(self.NUM_CASES_TAG,
                                                                                     str(len(case_set)))
            case_list_file.write("case_list_description: " + case_list_description + "\n")
            case_list_file.write("case_list_category: " + case_list_config_fields[
                case_list_config_header.index("META_CASE_LIST_CATEGORY")] + "\n")
            case_list_file.write("case_list_ids: " + '\t'.join(case_set) + "\n")

    def main(self):
        if self.verbose:
            print("LOG: case_list_config_file='%s'" % (self.config_file_name))
            print("LOG: case_list_dir='%s'" % (self.case_list_dir))
            print("LOG: study_dir='%s'" % (self.study_dir))
            print("LOG: study_id='%s'" % (self.study_id))
            print("LOG: overwrite='%s'" % (self.overwrite))
            print("LOG: verbose='%s'" % (self.verbose))

        if not isfile(self.config_file_name):
            print("ERROR: case list configuration file '%s' does not exist or is not a file" % (
                self.config_file_name),
                  file=sys.stderr)
            sys.exit(2)

        if not isdir(self.case_list_dir):
            print("ERROR: case list file directory '%s' does not exist or is not a directory" % (self.case_list_dir),
                  file=sys.stderr)
            sys.exit(2)

        if not isdir(self.study_dir):
            print("ERROR: study directory '%s' does not exist or is not a directory" % (self.study_dir),
                  file=sys.stderr)
            sys.exit(2)

        self.generate_case_lists()
