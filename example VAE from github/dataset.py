import re
from torchtext import data


class MyDataset(data.Dataset):

    def __init__(self, path, sequence_field, **kwargs):
        fields = [('sequence', sequence_field)]
        examples = []
        with open(path, 'r') as f:
            for sequence in f:
                examples.append(data.Example.fromlist([sequence], fields))
        super(MyDataset, self).__init__(examples, fields, **kwargs)

    @classmethod
    def splits(cls, sequence_field, train='train', **kwargs):

        return super(MyDataset, cls).splits(sequence_field=sequence_field, train=train, **kwargs)

def get_iterators(opt):
    sequence_field = data.Field(init_token=opt.start_token, eos_token=opt.end_token, lower=False, batch_first=True)
    train_data, val_data = MyDataset.splits(path="../data/spike_proteins", train ="1_in_500_cleaned.txt", test="1_in_500_cleaned.txt", sequence_field=sequence_field)
    sequence_field.build_vocab(train_data, val_data, max_size=opt.n_vocab-4, vectors='glove.6B.300d')
    train_vocab = sequence_field.vocab

    train_iter, val_iter = data.BucketIterator.splits((train_data, val_data), batch_size=opt.batch_size, sort_key = lambda x: len(x.text), repeat = False)
    return train_iter, val_iter, train_vocab