
from worker_genome import *

class AgeData:
    class INFO:
        def __init__(self, gene_id, gene_age, group_label):
            self.gene_id = gene_id
            self.gene_age = gene_age
            self.group_label = group_label

        def get_gene_age_numeric(self):
            if self.gene_age.__contains__(">"):
                return int(self.gene_age.replace(">", ""))
            return int(self.gene_age)

    MIN_SIZE_GROUP = 150

    def label_to_numeric(self, x):
        return int(x.replace(">", "")) if x.__contains__(">") else int(
            x.split("-")[0]) if x.__contains__("-") else int(x)

    def relabelling_gene_ages(self):
        while True:
            label_freqs = {}
            for gene_id, age_data in self.data_by_gene.items():
                if not label_freqs.__contains__(age_data.group_label):
                    label_freqs[age_data.group_label] = 0
                label_freqs[age_data.group_label] += 1

            small_group_label = ""
            arr = []
            for label, size in label_freqs.items():
                arr.append(label)
                if size < self.MIN_SIZE_GROUP:
                    small_group_label = label

            if small_group_label == "": break

            arr.sort(key=lambda x: self.label_to_numeric(x))

            partner_group_label = ""
            for i in range(0, len(arr)):
                if arr[i] == small_group_label:
                    if i == 0:
                        partner_group_label = arr[i + 1]
                    elif arr[i + 1].__contains__(">"):
                        partner_group_label = arr[i - 1]
                    else:
                        partner_group_label = arr[i + 1] if label_freqs[arr[i + 1]] < label_freqs[arr[i - 1]] else \
                            arr[i - 1]
                    break
            if self.label_to_numeric(small_group_label) < self.label_to_numeric(partner_group_label):
                new_label = f"{small_group_label.split('-')[0]}-{partner_group_label.split('-')[1]}"
            else:
                new_label = f"{partner_group_label.split('-')[0]}-{small_group_label.split('-')[1]}"

            for gene_id, age_data in self.data_by_gene.items():
                if age_data.group_label == small_group_label or age_data.group_label == partner_group_label:
                    age_data.group_label = new_label
                self.data_by_gene[gene_id] = age_data

    def __init__(self, spec: SPECIES, gen: GenomeWorker):
        self.spec = spec
        self.gen = gen

        self.data_by_gene = {}

        self.load_age_data()
        self.relabelling_gene_ages()

    def load_age_data(self):
        # http://genorigin.chenzxlab.cn/
        ages_file = open(f"./used_data/age_data/{self.spec.name}.csv", "r")
        lines = ages_file.readlines()

        for age_line in lines:
            if age_line.__contains__("ensembl_id"): continue
            arr = age_line.replace("\n", '').split(",")
            gene_id, gene_age, age_interval = arr[0], arr[1], arr[2]

            # we don't need gene age data, if gene is not loaded in genome
            if self.gen.feature_by_id(f"gene:{gene_id}") is None: continue

            data = self.INFO(gene_id, gene_age=gene_age, group_label=age_interval)
            assert not self.data_by_gene.__contains__(gene_id)
            self.data_by_gene[gene_id] = data

        print(f"Loaded AgeData (http://genorigin.chenzxlab.cn/) for {self.spec.short_name()}:\n\t"
              f"{len(self.data_by_gene.keys())} genes/records")

    class AgeGroupsType(enum.Enum):
        same_size_age_steps = 0,
        same_size_groups = 1,
        annotated_groups = 2

    gene_group_size = 800
    age_step_size = 50

    def calculate_groups(self, groups_type: AgeGroupsType):
        data = list(self.data_by_gene.values())
        groups = []

        if groups_type == self.AgeGroupsType.annotated_groups:
            data_by_annotated_label = {}
            for age_info in data:
                annotated_age_group_label = age_info.group_label
                if not data_by_annotated_label.__contains__(annotated_age_group_label):
                    data_by_annotated_label[annotated_age_group_label] = []
                data_by_annotated_label[annotated_age_group_label].append(age_info)
            arr = list(data_by_annotated_label.keys())
            arr.sort(key=lambda x: self.label_to_numeric(x))
            for interval in arr:
                groups.append(data_by_annotated_label[interval])

        if groups_type == self.AgeGroupsType.same_size_groups:
            data.sort(key=lambda x: x.get_gene_age_numeric())
            groups_count = (len(data) - 1) // self.gene_group_size + 1
            for index in range(groups_count):
                arr = []
                for i in range(index * self.gene_group_size, (index + 1) * self.gene_group_size):
                    if len(data) > i: arr.append(data[i])
                groups.append(arr)

        if groups_type == self.AgeGroupsType.same_size_age_steps:
            data.sort(key=lambda x: x.get_gene_age_numeric())
            for age_info in data:
                age = age_info.get_gene_age_numeric()
                index = (age - 1) // self.age_step_size
                while len(groups) <= index:
                    groups.append([])
                groups[index].append(age_info)
        return groups

    def get_group_label(self, group, group_index, age_group_type: AgeGroupsType):
        if age_group_type == self.AgeGroupsType.annotated_groups:
            assert len(group) > 0
            return group[0].group_label
        elif age_group_type == self.AgeGroupsType.same_size_groups:
            min_age, max_age = -1, -1
            for age_info in group:
                age = age_info.get_gene_age_numeric()
                min_age = age if min_age == -1 else min(age, min_age)
                max_age = age if max_age == -1 else max(age, max_age)
            return f"{min_age}-{max_age}"
        else:
            return f"{group_index * self.age_step_size + 1}-{(group_index + 1) * self.age_step_size}"

    def get_gene_ids_in_group(self, group):
        arr = []
        for age_info in group:
            gene_id = age_info.gene_id
            gene = self.gen.feature_by_id(f"gene:{gene_id}")
            assert gene is not None
            arr.append(f"gene:{gene_id}")
        return arr

    def get_gene_age(self, gene_id):
        gene_id = gene_id.replace("gene:", "")
        if self.data_by_gene.__contains__(gene_id):
            return self.data_by_gene[gene_id].get_gene_age_numeric()
        return None

