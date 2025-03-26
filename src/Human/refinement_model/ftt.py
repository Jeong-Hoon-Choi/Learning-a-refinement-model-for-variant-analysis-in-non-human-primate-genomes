import torch
from tab_transformer_pytorch import FTTransformer
from torch.utils.data import DataLoader, TensorDataset
import torch.optim as optim
from sklearn.model_selection import train_test_split
import torch.nn.functional as F
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score
import time


def FTT(train_feature, train_label, test_feature, test_label):
    model = FTTransformer(
        categories=[],
        num_continuous=train_feature.shape[1],
        dim=128,
        dim_out=2,
        depth=6,
        heads=8,
        attn_dropout=0.1,
        ff_dropout=0.1
    )

    device = ('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(device)

    model = model.to(device)

    train_features, val_features, train_labels, val_labels = train_test_split(train_feature, train_label, test_size=0.2, random_state=42)

    train_dataset = TensorDataset(torch.tensor(train_features.values, dtype=torch.float32), torch.tensor(train_labels.values, dtype=torch.int64))
    val_dataset = TensorDataset(torch.tensor(val_features.values, dtype=torch.float32), torch.tensor(val_labels.values, dtype=torch.int64))
    test_dataset = TensorDataset(torch.tensor(test_feature.values, dtype=torch.float32), torch.tensor(test_label.values, dtype=torch.int64))

    train_loader = DataLoader(train_dataset, batch_size=256, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=256, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=256, shuffle=False)

    optimizer = optim.Adam(model.parameters(), lr=0.001)
    criterion = torch.nn.CrossEntropyLoss()

    epochs = 10
    for epoch in range(epochs):
        time_i = time.time()
        model.train()
        for batch_data, batch_labels in train_loader:
            optimizer.zero_grad()
            empty_categories = torch.empty((batch_data.shape[0], 0)).to(device)
            outputs = model(empty_categories, batch_data.to(device))
            loss = criterion(outputs, batch_labels.to(device))
            loss.backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            val_loss = 0.0
            for batch_data, batch_labels in val_loader:
                empty_categories = torch.empty((batch_data.shape[0], 0)).to(device)
                outputs = model(empty_categories, batch_data.to(device))
                val_loss += criterion(outputs, batch_labels.to(device)).item()
            val_loss /= len(val_loader)

        time_now = round(time.time() - time_i, 2)
        print(f"Epoch {epoch + 1}/{epochs}, Validation Loss: {val_loss:.4f}", 'time :', time_now)

    model.eval()

    all_probabilities = []
    all_predictions = []

    with torch.no_grad():
        for batch_data, batch_labels in test_loader:
            empty_categories = torch.empty((batch_data.shape[0], 0)).to(device)
            outputs = model(empty_categories, batch_data.to(device))
            probabilities = F.softmax(outputs, dim=1)  # 확률을 계산
            all_probabilities.extend(probabilities.cpu().numpy())

            _, predicted = outputs.max(1)
            all_predictions.extend(predicted.cpu().numpy())

    return model, all_predictions, all_probabilities


def FFT_inf(model_dir, test_feature, test_label):
    model = FTTransformer(
        categories=[],
        num_continuous=test_feature.shape[1],
        dim=128,
        dim_out=2,
        depth=6,
        heads=8,
        attn_dropout=0.1,
        ff_dropout=0.1
    )

    device = ('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)

    model.load_state_dict(torch.load(model_dir))
    model.eval()

    test_dataset = TensorDataset(torch.tensor(test_feature.values, dtype=torch.float32), torch.tensor(test_label.values, dtype=torch.int64))
    test_loader = DataLoader(test_dataset, batch_size=256, shuffle=False)

    all_probabilities = []
    all_predictions = []

    with torch.no_grad():
        for batch_data, batch_labels in test_loader:
            empty_categories = torch.empty((batch_data.shape[0], 0)).to(device)
            outputs = model(empty_categories, batch_data)
            probabilities = F.softmax(outputs, dim=1)  # 확률을 계산
            all_probabilities.extend(probabilities.cpu().numpy())

            _, predicted = outputs.max(1)
            all_predictions.extend(predicted.cpu().numpy())

    return all_probabilities, all_predictions
