import os
import sys

class PosGroups:
    @staticmethod
    def group_pos(dicts):
        def similar_lists(list1,list2):
            return list1 == list2
        def get_lists(dicts):
            lists = []
            sample_counts = {}
            for type_val in dicts.values():
                for group_val in type_val.values():
                    for condition,condition_val in group_val.items():
                        for sample,sample_list in condition_val.items():
                            lists.append(sample_list)
                            list_tuple = tuple(sample_list)
                            if list_tuple not in sample_counts:
                                sample_counts[list_tuple] = {}
                            if condition not in sample_counts[list_tuple]:
                                sample_counts[list_tuple][condition] = set()
                            sample_counts[list_tuple][condition].add(sample)

            return lists,sample_counts

        def group_lists(lists):
            groups = []
            for lst in lists:
                matched = False
                for group in groups:
                    if lst == group[0]:
                        group.append(lst)
                        matched = True
                        break
                if not matched:
                    groups.append([lst])
            return groups

        all_lists,sample_counts = get_lists(dicts)
        grouped_lists = group_lists(all_lists)
        num_groups = len(grouped_lists)
        unique_lists = [(group[0], {condition: len(samples) for condition, samples in sample_counts[tuple(group[0])].items()}) for group in grouped_lists]

        return num_groups, unique_lists

