#!/usr/bin/env python
"""What is this?."""

import struct


def parseIPCMessage(data):
    """TODO DOCUMENTME."""
    (bodyLength, nameLength) = struct.unpack("<HB", data[:3])

    if not bodyLength:
        return

    name = data[3 : 3 + nameLength]
    valType = struct.unpack("<B", data[3 + nameLength : 3 + nameLength + 1])[0]
    value = 0

    if valType == 0:  # long
        value = struct.unpack(
            "<q", data[3 + nameLength + 1 : 3 + nameLength + 1 + 8]
        )[0]
    elif valType == 1:  # double
        value = struct.unpack(
            "<d", data[3 + nameLength + 1 : 3 + nameLength + 1 + 8]
        )[0]

    return name, value


# class IPCThreadedTCPRequestHandler(socketserver.BaseRequestHandler):
# def __init__(self):
# SocketServer.BaseRequestHandler.__init__(self)

# def setup(self):
# print('ipcserver: got new connection %d from %s' %
# (self.request.fileno(), self.client_address))
# welcome = "Welcome on ipc-server. Client " + \
# str(self.request.fileno()) #+ str(self.server_address)
# out_msg = struct.pack("<B", len(welcome)) + \
# welcome + struct.pack("<B", 0) + struct.pack("<q", -111)
# self.request.send(struct.pack("<H", len(out_msg)) + out_msg)

# def handle(self):
# while 1:
# data = self.request.recv(1024)

# if data:
# cur_thread = threading.currentThread()
# start = 0
# end = 0
# while start < len(data):
# end = start + struct.unpack("<H",
# data[ start : start + 2 ])[0] +2
# print "parsing:", len(data), start, end
# name,value = parseIPCMessage(data[ start : end ])
# print("server: client(", self.request.fileno(), ")",
# name, value, "\n")
# start = end
# just echoing
# self.request.send(data)
# continue
# else:
# break;
# self.request.send(response)

# def finish(self):
# print("unsubscribe: ", self.request.fileno())

# class IPCServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
# pass
