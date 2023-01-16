import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from torch.optim.lr_scheduler import ReduceLROnPlateau
import random
from scipy.stats import spearmanr
import os

random.seed(42)
torch.manual_seed(42)
np.random.seed(42)

class AbacusNet(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.layer0 = nn.Sequential(
            nn.Linear(5,128),
            nn.GELU(),
            nn.Linear(128, 32),
            nn.GELU(),
            nn.Linear(32, 1)
            )
    def forward(self, x):
        return self.layer0(x)

class myDataset(Dataset):
    def __init__(self, abacus_df) -> None:
        super().__init__()
        self.data_df = abacus_df[['sai', 's1', 's2','pack', 'hb', 'ddG']]

    def __len__(self):
        return len(self.data_df)

    def __getitem__(self, index):
        return self.data_df.iloc[index].values[:5].astype(np.float32), self.data_df.iloc[index].values[-1].astype(np.float32)

def setup_train(data_csv='../train_abacus2.csv', cv_num=6):
    df = pd.read_csv(data_csv)
    batch_size = 128
    net = AbacusNet()
    cal_loss = nn.MSELoss()
    opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
    scheduler = ReduceLROnPlateau(opt, 'min')
    train_df = df[df['group'] != cv_num]
    test_df = df[df['group'] == cv_num]
    train_data = myDataset(train_df)
    test_data = myDataset(test_df)
    train_loader = DataLoader(train_data, shuffle=True, batch_size=batch_size)
    test_loader = DataLoader(test_data, shuffle=False, batch_size=512)
    return net, opt, train_loader, test_loader, cal_loss, scheduler

def run_nn_train(net, opt, train_loader, test_loader, cal_loss, scheduler):
    res_dict = {
        "train_losses" : [],
        "train_r" : [],
        "test_losses" : [],
        "test_r" : [],
        "lrs" : []
    }
    train_losses = []
    train_r = []
    test_losses = []
    test_r = []
    lrs = []
    for epoch in range(50):
        train_pred = []
        train_y = []
        train_loss_all = 0
        for x, y in train_loader:
            opt.zero_grad()
            pred = net(x)
            train_loss = cal_loss(pred.ravel(), y)
            train_loss.backward()
            opt.step()
            train_pred += pred.detach().ravel().tolist()
            train_y += y.ravel().tolist()
            train_loss_all += train_loss.item() * len(x)
        train_losses.append(train_loss_all/len(train_loader.dataset))
        with torch.no_grad():
            test_pred = []
            test_y = []
            for x, y in test_loader:
                pred = net(x)
                test_loss = cal_loss(pred.ravel(), y)
                test_pred += pred.ravel().tolist()
                test_y += y.ravel().tolist()
                test_losses.append(test_loss.item())
            scheduler.step(test_loss)
            lr = scheduler._last_lr[0]
            lrs.append(lr)
        # if (epoch + 1)%10 == 0:
        #     print(f"Epoch {epoch+1}")
        #     print("Train r: {:.4f}, Test r: {:.4f}".format(spearmanr(train_pred, train_y)[0], spearmanr(test_pred, test_y)[0]))
        train_r.append(spearmanr(train_pred, train_y)[0])
        test_r.append(spearmanr(test_pred, test_y)[0])
    res_dict['train_losses'] = train_losses
    res_dict["train_r"] = train_r
    res_dict["test_losses"] = test_losses
    res_dict["test_r"] = test_r
    res_dict["lrs"] = lrs
    return res_dict, net, test_pred, test_y

def get_models(path_to_checkpoints=''):
    models = []
    for i in range(10):
        net = AbacusNet()
        net.load_state_dict(torch.load(os.path.join(os.path.expanduser('~'),f".cache/ddgscan/abacus2_nn_{i}.pt")))
        models.append(net)
    return models

# model_list = get_models('/home/jsun/egnn-ddg-predictor/DDGScan_DL/abacus2_data')
# with torch.no_grad():
#     all_preds = []
#     for net in model_list:
#         x = torch.tensor(fp_abacus_df[['sai', 's1', 's2', 'pack', 'hb']].values).float()
#         pred0 = net(x)
#         pred0 = pred0.ravel().numpy()
#         all_preds.append(pred0)
